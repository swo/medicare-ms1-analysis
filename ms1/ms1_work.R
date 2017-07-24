# Do all the heavy-lifting analyses so that the Rmd can just show the results

sum_over = function(target_columns) {
  dummy_vars = sapply(1:length(target_columns), function(i) sprintf('x%i', i))
  dummy_vars_string = paste0(dummy_vars, collapse=',')
  dummy_vars_list = lapply(target_columns, as.name) %>% setNames(dummy_vars)
  f = as.formula(paste0(c('~sum(', dummy_vars_string, ')'), collapse=''))
  lazyeval::interp(f, .values=dummy_vars_list)
}

# Load the census regions. Code them as factors so that Northeast is taken as
# the baseline in the linear models.
# NB: I put DC into the South
regions = read_tsv('../../db/census-regions/census-regions.tsv') %>%
  select(state, region)

# Keep chronic conditions that have a single reference year. Exclude Alzheimer's,
# which is a subset of Alz. & dementia
condition_names = read_tsv('../chronic-conditions/names.tsv') %>% filter(condition_code != 'ALZH')
names_1yr_ccs = condition_names %>% filter(condition_ref_years==1) %$% condition_code

dx_codes = read_tsv('../../data/fd_codes.tsv') %>%
  select(code, diagnosis_type)

dx_types = unique(dx_codes$diagnosis_type)

# window of days that you can look ahead from an rx to a dx
rxdx_window = 3

load_data = function(year) {
  # pde data
  pde = read_tsv(sprintf('../pde_%i.tsv', year)) %>%
    select(bene_id, service_date, antibiotic) %>%
    distinct() %>%
    mutate(year=year, pde_id=1:n())

  # diagnosis claims
  dx = read_tsv(sprintf('../dx_%i.tsv', year)) %>%
    mutate(year=year)
  
  # associated each PDE with a diagnosis, if possible
  pde_dx = pde %>%
    left_join(dx, by=c('year', 'bene_id')) %>%
    filter(between(service_date - from_date, 0, rxdx_window)) %>%
    group_by(pde_id) %>% filter(from_date==max(from_date)) %>% ungroup() %>%
    select(pde_id, diagnosis_type) %>%
    right_join(pde, by='pde_id') %>%
    replace_na(list(diagnosis_type='none'))

  # total numbers of PDEs
  n_claims = count(pde, bene_id) %>% rename(n_claims=n)

  # count chronic conditions. 'n_cc' means one-year cc's
  cc = read_feather(sprintf('../cc_%i.feather', year)) %>%
    rowwise() %>%
    mutate_(n_cc=sum_over(names_1yr_ccs)) %>%
    ungroup()

  # join the cc, summary PDE, and dx data into bene
  bene = read_tsv(sprintf('../bene_%i.tsv', year)) %>%
    mutate(is_female=sex=='female', is_white=race=='white', is_dual=buyin_months>0) %>%
    left_join(regions, by='state') %>%
    mutate(year=year) %>%
    left_join(cc, by='bene_id') %>%
    left_join(n_claims, by='bene_id') %>% replace_na(list(n_claims=0))

  # take the bene, pde (with diagnoses) and diagnoses (separately)
  list(bene=bene, pde=pde_dx, dx=dx)
}

dat = lapply(2011:2014, load_data)

bene = lapply(dat, function(df) df$bene) %>% bind_rows %>%
  group_by(bene_id) %>%
  mutate(n_years=n(), in_cohort=n_years==4) %>%
  ungroup() %>%
  filter(in_cohort) %>% select(-in_cohort)

cohort_ids = unique(bene$bene_id)

pde = lapply(dat, function(df) df$pde) %>% bind_rows %>%
  filter(bene_id %in% cohort_ids)

dx = lapply(dat, function(df) df$dx) %>% bind_rows() %>%
  filter(bene_id %in% cohort_ids)

output_table = function(x, base) write_tsv(x, sprintf('tables/%s.tsv', base))
sem = function(x) sd(x) / sqrt(length(x))

summarize_totals = function(df) {
  summarize(df, n_bene=n(),
            mean_age=mean(age),
            sd_age=sd(age),
            sem_age=sem(age),
            mean_n_cc=mean(n_cc),
            sd_n_cc=sd(n_cc),
            sem_n_cc=sem(n_cc),
            n_female=sum(is_female),
            n_white=sum(is_white),
            n_dual=sum(is_dual),
            n_region_south=sum(region=='South'),
            n_region_west=sum(region=='West'),
            n_region_northeast=sum(region=='Northeast'),
            n_region_midwest=sum(region=='Midwest'))
}

# summary characteristics
totals = bene %>%
  group_by(year) %>%
  summarize_totals() %T>%
  output_table('totals')

# the top few inidividual antibiotics
claims_by_abx = pde %>%
  count(year, antibiotic) %>% ungroup() %>%
  left_join(select(totals, year, n_bene), by='year') %>%
  mutate(cpkp=n*1000/n_bene) %T>%
  output_table('claims_by_abx')

# keep track of which are the top 10 abx
top_abx = pde %>%
  count(antibiotic) %>%
  arrange(desc(n)) %$%
  head(antibiotic, 10)

# from here on, the denominators are taken from the beneficiary data grouped
# in the same way as the consumption data, so we can use the summarize_by function
summarize_by = function(bene, by) {
  group_by_(bene, 'year', by) %>%
    summarize(n_bene=n(), n_claims=sum(n_claims)) %>%
    ungroup() %>%
    mutate(cpkp=n_claims*1000/n_bene) %>%
    output_table(paste0('claims_by_', by))
}

summarize_by(bene, 'sex')
summarize_by(bene, 'region')
summarize_by(bene, 'race')

# then by age
# (the input to summarize_by is modified, since we want to group by age *group*,
# not just by raw age)
bene %>%
  mutate(age_group=case_when(
    between(.$age, 66, 76) ~ '66-76',
    between(.$age, 77, 86) ~ '77-86',
    between(.$age, 87, 96) ~ '87-96')) %>%
  summarize_by('age_group')

# run a model, tidy, and save the output
subtract_min = function(x) x - min(x)
model_f = function(df, frmla) {
  dat = df %>% mutate_at(vars(year, age), subtract_min)
  f = function(...) glm(formula=frmla, data=dat, ...)
  model_lm = f()
  model_poisson = f(family='poisson')
  model_log = f(family=gaussian(link='log'), start=model_poisson$coefficients)
  
  bind_rows(
    mutate(tidy(model_lm), model='lm'),
    mutate(tidy(model_poisson), model='poisson'),
    mutate(tidy(model_log), model='log')
  )
}

# models predicting consumption
# models with all beneficiaries
model_f(bene, n_claims ~ year) %>% output_table('model_total_unadjusted')
model_f(bene, n_claims ~ year + age + n_cc + is_dual + region) %>% output_table('model_total_adjusted')

# models for each individual drug
single_abx_model = function(abx, frmla) {
  pde %>%
    filter(antibiotic==abx) %>%
    count(year, bene_id) %>% ungroup() %>% rename(y=n) %>%
    right_join(bene, by=c('year', 'bene_id')) %>%
    replace_na(list(y=0)) %>%
    model_f(frmla) %>%
    mutate(antibiotic=abx)
}

model_abx_overall = bene %>%
  rename(y=n_claims) %>%
  model_f(y ~ year + age + n_cc + is_dual + region) %>%
  mutate(antibiotic='overall')

model_abx = lapply(top_abx, function(a) single_abx_model(a, y ~ year + age + n_cc + is_dual + region)) %>%
  bind_rows() %>%
  bind_rows(model_abx_overall) %T>%
  output_table('model_abx')

# diagnoses analysis
times = data_frame(year=c(2011, 2014), time=c(0, 1))
time_f = function(df) left_join(times, df, by='year')

## age structure
p_age = bene %>% time_f() %>%
  count(time, age) %>%
  group_by(time) %>% mutate(p_age=n/sum(n)) %>% ungroup() %>%
  mutate(key=sprintf('p_age_t%i', time)) %>%
  select(age, key, p_age) %>%
  spread(key, p_age)

## a "narrow" age structure
# it excludes the old and young ages that don't appear in both
# the 2011 and 2014 data, so we can adjust for age structure
p_age_narrow = p_age %>%
  select(age, p_age=p_age_t0) %>%
  filter(between(age, 69, 92)) %>%
  mutate(p_age=p_age/sum(p_age))

n_dxage = dx %>%
  left_join(select(bene, year, bene_id, age), by=c('year', 'bene_id')) %>%
  time_f() %>%
  count(time, age, diagnosis_type) %>% ungroup()

n_rxdxage = pde %>%
  left_join(select(bene, year, bene_id, age), by=c('year', 'bene_id')) %>%
  time_f() %>%
  count(time, age, antibiotic, diagnosis_type) %>% ungroup()

f_rxdx = n_rxdxage %>%
  group_by(antibiotic, diagnosis_type) %>%
  summarize(n=sum(n)) %>%
  group_by(antibiotic) %>%
  mutate(f_rxdx=n/sum(n)) %>%
  ungroup() %>%
  select(abx=antibiotic, diagnosis_type, f_rxdx)

proportion_f = function(abx) {
  dx_counts = n_dxage %>%
    mutate(key=sprintf('n_dx_t%i', time)) %>%
    select(age, diagnosis_type, key, n) %>%
    spread(key, n)
  rx_counts = n_rxdxage %>%
    filter(antibiotic==abx) %>%
    mutate(key=sprintf('n_rx_t%i', time)) %>%
    select(age, diagnosis_type, key, n) %>%
    spread(key, n)
  
  overall_rx_counts = n_rxdxage %>%
    filter(antibiotic==abx) %>%
    group_by(time, age) %>% summarize(n=sum(n)) %>% ungroup() %>%
    mutate(key=sprintf('n_rx_t%i', time)) %>%
    select(age, key, n) %>%
    spread(key, n) %>% mutate(r=n_rx_t1/n_rx_t0) %>%
    left_join(p_age, by='age') %>%
    mutate(r=(n_rx_t1/n_rx_t0)/(p_age_t1/p_age_t0)) %>%
    select(age, r) %>%
    mutate(diagnosis_type='overall')
  
  inner_join(dx_counts, rx_counts, by=c('age', 'diagnosis_type')) %>%
    mutate(r=(n_rx_t1/n_rx_t0)/(n_dx_t1/n_dx_t0)) %>%
    select(age, diagnosis_type, r) %>%
    bind_rows(overall_rx_counts) %>%
    filter(between(age, 69, 92)) %>%
    left_join(p_age_narrow, by='age') %>%
    group_by(diagnosis_type) %>%
    summarize(r=sum(r*p_age))
}

rxdx_result = lapply(top_abx, function(x) {
  proportion_f(x) %>%
    select(diagnosis_type, value=r) %>%
    mutate(abx=x)
}) %>% bind_rows() %>%
  left_join(f_rxdx, by=c('abx', 'diagnosis_type')) %>%
  mutate(f_rxdx=if_else(diagnosis_type=='overall', 1.0, f_rxdx)) %>%
  select(abx, diagnosis_type, p_rxdx_ratio=value, f_rxdx) %>%
  output_table('rxdx')
