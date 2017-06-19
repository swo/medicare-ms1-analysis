# Do all the heavy-lifting analyses so that the Rmd can just show the results

import::from(broom, tidy)

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

condition_names = read_tsv('../chronic-conditions/names.tsv') %>% filter(condition_code != 'ALZH')
names_1yr_ccs = condition_names %>% filter(condition_ref_years==1) %$% condition_code

dx_codes = read_tsv('../sex-difference/sex_codes.tsv') %>%
  select(code, diagnosis_type)

load_data = function(year) {
  # pde data
  pde = read_tsv(sprintf('../pde_%i.tsv', year)) %>%
    select(bene_id, antibiotic, days_supply) %>%
    group_by(bene_id, antibiotic) %>%
    summarize(n_claims=n(), days_supply=sum(days_supply)) %>%
    ungroup() %>%
    mutate(year=year)
  
  # diagnosis claims
  dx_car_fn = sprintf('../../data/dx_car_claims_%i.tsv', year)
  dx_op_fn = sprintf('../../data/dx_op_claims_%i.tsv', year)
  
  dx = bind_rows(read_tsv(dx_car_fn),
                 read_tsv(dx_op_fn)) %>%
    mutate(date=dmy(from_date)) %>%
    rename(bene_id=BENE_ID, code=diagnosis) %>%
    left_join(dx_codes, by='code') %>%
    filter(diagnosis_type=='acute_rc') %>%
    select(bene_id, date) %>%
    group_by(bene_id) %>%
    summarize(n_acute_rc_dx=length(unique(date)))

  # total numbers of PDEs
  n_claims = pde %>% group_by(bene_id) %>% summarize(n_claims=sum(n_claims))

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
    left_join(n_claims, by='bene_id') %>% replace_na(list(n_claims=0)) %>%
    left_join(dx, by='bene_id') %>% replace_na(list(n_acute_rc_dx=0))

  list(bene=bene, pde=pde)
}

dat = lapply(2011:2014, load_data)

bene = lapply(dat, function(df) df$bene) %>% bind_rows %>%
  group_by(bene_id) %>%
  mutate(n_years=n(),
         start_age=min(age),
         in_cohort=n_years==4) %>%
  ungroup()

pde = lapply(dat, function(df) df$pde) %>% bind_rows %>%
  left_join(distinct(select(bene, bene_id, in_cohort)), by='bene_id')

#rm(dat)

output_table = function(x, base) write_tsv(x, sprintf('tables/tbl_%s.tsv', base))
sem = function(x) sd(x) / length(x)

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

cohort_totals = bene %>%
  filter(in_cohort) %>%
  group_by(year) %>%
  summarize_totals() %T>%
  output_table('cohort_totals')

# the top few inidividual antibiotics
claims_by_abx = pde %>%
  group_by(antibiotic, year) %>%
  summarize_at(vars(n_claims, days_supply), sum) %>%
  ungroup() %>%
  left_join(select(totals, year, n_bene), by='year') %>%
  mutate(cpkp=n_claims*1000/n_bene, did=days_supply*1000/(n_bene*365)) %T>%
  output_table('claims_by_abx')

claims_by_abx = pde %>%
  filter(in_cohort) %>%
  group_by(antibiotic, year) %>%
  summarize_at(vars(n_claims, days_supply), sum) %>%
  ungroup() %>%
  left_join(select(cohort_totals, year, n_bene), by='year') %>%
  mutate(cpkp=n_claims*1000/n_bene, did=days_supply*1000/(n_bene*365)) %T>%
  output_table('cohort_claims_by_abx')

# keep track of which are the top 10 abx
top_abx = claims_by_abx %>%
  group_by(antibiotic) %>%
  summarize(cpkp=mean(cpkp)) %>%
  arrange(desc(cpkp)) %>%
  head(10) %$%
  antibiotic

# from here on, the denominators are taken from the beneficiary data grouped
# in the same way as the consumption data, so we can use the summarize_by function
summarize_by = function(bene, by) {
  group_by_(bene, 'year', by) %>%
    summarize(n_bene=n(), n_claims=sum(n_claims)) %>%
    ungroup() %>%
    mutate(cpkp=n_claims*1000/n_bene)
}

summarize_by_with_cohort = function(bene, by, table_name) {
  summarize_by(bene, by) %>% output_table(table_name)
  bene %>% filter(in_cohort) %>% summarize_by(by) %>% output_table(paste0('cohort_', table_name))
}

summarize_by_with_cohort(bene, 'sex', 'claims_by_sex')
summarize_by_with_cohort(bene, 'region', 'claims_by_region')
summarize_by_with_cohort(bene, 'race', 'claims_by_race')

# then by age
# (the input to summarize_by is modified, since we want to group by age *group*,
# not just by raw age)
bene %>%
  mutate(age_group=case_when(
    between(.$age, 66, 76) ~ '66-76',
    between(.$age, 77, 86) ~ '77-86',
    between(.$age, 87, 96) ~ '87-96')) %>%
  summarize_by_with_cohort('age_group', 'claims_by_age')

# run a model, tidy, and save the output
model_f = function(bene, frmla, table_name) {
  lm(formula=frmla, data=bene) %>% tidy %>% output_table(table_name)
}

# models predicting consumption
# models with all beneficiaries
model_f(bene, n_claims ~ year, 'model1')
model_f(bene, n_claims ~ year + age*n_cc + is_female + is_dual + is_white + region, 'model2')

# models using just the cohort
lm_cohort_bene = bene %>% filter(in_cohort)

model_f(lm_cohort_bene, n_claims ~ year, 'cohort_model1')
model_f(lm_cohort_bene, n_claims ~ year + age*n_cc + is_female + is_dual + is_white + region, 'cohort_model2')

# models for each individual drug
subtract_min = function(x) x - min(x)
single_abx_model = function(bene, abx, frmla, y_name='n_claims') {
  pde %>%
    rename_(y=y_name) %>% 
    filter(antibiotic==abx) %>%
    right_join(bene, by=c('year', 'bene_id')) %>%
    replace_na(list(y=0)) %>%
    mutate_at(vars(year, age), subtract_min) %>%
    lm(formula=frmla, data=.) %>%
    tidy %>%
    mutate(model=abx)
}

f = function(a) single_abx_model(bene, a, y ~ year + age*n_cc + is_female + is_dual + is_white + region)
lapply(top_abx, f) %>%
  bind_rows() %T>%
  output_table('model_abx')

f = function(a) single_abx_model(lm_cohort_bene, a, y ~ year + age*n_cc + is_female + is_dual + is_white + region)
lapply(top_abx, f) %>%
  bind_rows() %T>%
  output_table('cohort_model_abx')

# special levo & azithro models
azithro_f = y ~ year + age + is_female + is_dual + is_white + region + COPD + COPD:year + COPD:age + n_acute_rc_dx + n_acute_rc_dx:year + n_acute_rc_dx:age
single_abx_model(lm_cohort_bene, 'azithromycin', azithro_f) %>%
  output_table('cohort_model_azithro')
levo_f = y ~ year + age + is_female + is_dual + is_white + region + n_acute_rc_dx + n_acute_rc_dx:year + n_acute_rc_dx:age
single_abx_model(lm_cohort_bene, 'levofloxacin', levo_f) %>%
  output_table('cohort_model_levo')

# repeated those models with days supply
single_abx_model(lm_cohort_bene, 'azithromycin', azithro_f, y_name='days_supply') %>%
  output_table('cohort_model_azithro_days')
single_abx_model(lm_cohort_bene, 'levofloxacin', levo_f) %>%
  output_table('cohort_model_levo_days')
