# Do all the heavy-lifting analyses so that the Rmd can just show the results
library(forcats)

# Load the census regions. Code them as factors so that Northeast is taken as
# the baseline in the linear models.
# NB: I put DC into the South
regions = read_tsv('../../db/census-regions/census-regions.tsv') %>%
  select(state, region)

# Keep chronic conditions that have a single reference year. Exclude Alzheimer's,
# which is a subset of Alz. & dementia
condition_names = read_tsv('../chronic-conditions/names.tsv') %>% filter(condition_code != 'ALZH')
names_1yr_ccs = condition_names %>% filter(condition_ref_years==1) %$% condition_code
sum_1yr_cc = paste0(names_1yr_ccs, collapse=' + ')

dx_codes = read_tsv('../../data/fd_codes.tsv') %>%
  select(code, diagnosis_type)

dx_types = unique(dx_codes$diagnosis_type)

# window of days that you can look ahead from an rx to a dx
rxdx_window = 7

load_data = function(year) {
  # pde data
  pde = read_tsv(sprintf('../pde_%i.tsv', year)) %>%
    select(bene_id, service_date, antibiotic) %>%
    distinct() %>%
    mutate(year=year, pde_id=1:n())

  # total numbers of PDEs
  n_claims = count(pde, bene_id) %>% rename(n_claims=n)

  # count chronic conditions. 'n_cc' means one-year cc's
  cc = read_feather(sprintf('../cc_%i.feather', year)) %>%
    mutate_(n_cc=sum_1yr_cc) %>%
    ungroup()

  # join the cc, summary PDE, and dx data into bene
  bene = read_tsv(sprintf('../bene_%i.tsv', year)) %>%
    mutate(year=year) %>%
    left_join(cc, by='bene_id') %>%
    left_join(n_claims, by='bene_id') %>% replace_na(list(n_claims=0))
  
  # diagnosis claims
  dx = read_tsv(sprintf('../../data/tmp_hcpcs_%i.tsv', year)) %>%
    select(bene_id=BENE_ID, code=dx, from_date=dt1) %>%
    # filter for bene's we have
    semi_join(bene, by='bene_id') %>%
    mutate(from_date=dmy(from_date), year=year, dx_id=1:n())
  
  # associate a PDE with each diagnosis, if possible
  dx_pde = dx %>%
    left_join(pde, by=c('year', 'bene_id')) %>%
    mutate(delay=service_date - from_date) %>%
    filter(between(delay, 0, rxdx_window)) %>%
    select(dx_id, antibiotic, delay) %>%
    right_join(dx, by='dx_id') %>%
    replace_na(list(antibiotic='none'))

  # take the bene, pde (with diagnoses) and diagnoses (separately)
  list(bene=bene, pde=pde, dx=dx_pde)
}

dat = lapply(2011:2014, load_data)

subtract_min = function(x) x - min(x)

bene = lapply(dat, function(df) df$bene) %>% bind_rows %>%
  mutate(heart_disease=AMI | ATRIALFB | CHF | ISCHMCHT) %>%
  mutate(is_female=sex=='female', is_white=race=='white', is_dual=buyin_months>0) %>%
  mutate(age=age-1) %>%
  left_join(regions, by='state') %>%
  mutate(region=factor(region, levels=c('Northeast', 'West', 'Midwest', 'South'))) %>%
  mutate(year0=subtract_min(year),
         age0=subtract_min(age))

pde = lapply(dat, function(df) df$pde) %>% bind_rows

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
summarize_by(bene, 'age')

# run a model and tidy
glm_f = function(df, frmla, ...) eval(substitute(function(df2) glm(frmla, data=df2, ...)))(df)
model_f = function(df, frmla, name) {
  glm_f(df, frmla, family='poisson') %>%
    tidy %>%
    mutate(name=name)
}

abx_frmla = y ~ year0 + age0 + n_cc + is_female + is_white + is_dual + region

# get models and predictions for all drugs combined
overall_model = bene %>%
  rename(y=n_claims) %>%
  model_f(abx_frmla, 'overall')

# get models for individual abx
single_f = function(abx, frmla) {
  pde %>%
    filter(antibiotic==abx) %>%
    count(year, bene_id) %>% ungroup() %>% rename(y=n) %>%
    right_join(bene, by=c('year', 'bene_id')) %>%
    replace_na(list(y=0)) %>%
    model_f(frmla, abx)
}

single_models = lapply(top_abx, function(a) single_f(a, abx_frmla))

# report the model parameters
bind_rows(overall_model, single_models) %>%
  output_table('model_abx')

# trends in all fluoroquinolones
model_fq = pde %>%
  filter(str_detect(antibiotic, 'floxacin')) %>%
  count(year, bene_id) %>% ungroup() %>% rename(y=n) %>%
  right_join(bene, by=c('year', 'bene_id')) %>% replace_na(list(y=0)) %>%
  model_f(abx_frmla, 'fluoroquinolone') %T>%
  output_table('model_fq')

# special heart disease models
# how many beneficiaries have heart disease?
bene %>%
  count(year, heart_disease) %>%
  output_table('counts_heartdisease')

# did azithro consumption decrease faster in benes with heart disease?
heart_frmla = y ~ year0*heart_disease + age0 + n_cc + is_female + is_white + is_dual + region
model_azithro = single_f('azithromycin', heart_frmla) %T>%
  output_table('model_azithro')

model_levo = single_f('levofloxacin', heart_frmla) %T>%
  output_table('model_levo')

## DIAGNOSES ##
dx_raw = lapply(dat, function(x) x$dx) %>% bind_rows() %>%
  left_join(dx_codes, by='code')

# count diagnoses
dx_counts = dx_raw %>%
  select(dx_id, year, diagnosis_type) %>%
  distinct() %>%
  count(year, diagnosis_type)

# count delays
dx_delays = dx_raw %>%
  count(delay) %T>%
  output_table('dx_delay')

# for each diagnosis, how often is each drug used?
dx_counts = dx_raw %>%
  count(year, diagnosis_type, antibiotic) %>% ungroup() %T>%
  output_table('dxrx_fractions')

# create a table of dx (rows) with multiple abx as columns
dxrx = dx_raw %>%
  select(year, bene_id, from_date, diagnosis_type, antibiotic) %>%
  mutate(antibiotic=fct_other(antibiotic, keep=c('none', top_abx))) %>%
  distinct() %>%
  mutate(dummy=TRUE) %>% spread(antibiotic, dummy, fill=FALSE) %>%
  left_join(bene, by=c('year', 'bene_id'))
  
dxrx_abx = c('azithromycin', 'levofloxacin', 'amoxicillin/clavulanate', 'cephalexin', 'ciprofloxacin')

# for a diagnosis type dx and an antibiotic, what is the trend in the probability
# that that drug will be used for that diagnosis?
# i.e., logistic regression drug_used? ~ dx_present?*year + covariates
# we also ask about "any diagnosis", i.e., drug_used? ~ year + covariates
covariates = str_c(c('age0', 'is_white', 'is_dual', 'is_female', 'region'), collapse=' + ')
dxrx_f = function(dxrx, dx, abx) {
  if (dx=='any_dx') {
    frmla = as.formula(str_interp("`${abx}` ~ year0 + ${covariates}"))
  } else {
    frmla = as.formula(str_interp("`${abx}` ~ dx*year0 + ${covariates}"))
  }
  
  dxrx %>%
    mutate(dx=diagnosis_type==dx) %>%
    glm_f(frmla, family='binomial') %>%
    tidy %>%
    mutate(diagnosis_type=dx, antibiotic=abx)
}


dxrx_res = crossing(diagnosis_type=c('any_dx', unique(dxrx$diagnosis_type)),
                    antibiotic=c(top_abx, 'none', 'Other')) %>%
  rowwise() %>%
  do(dxrx_f(dxrx, .$diagnosis_type, .$antibiotic)) %>%
  ungroup() %T>%
  output_table('model_dxrx')
