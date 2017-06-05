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

load_data = function(year) {
  bene = read_tsv(sprintf('../bene_%i.tsv', year)) %>%
    mutate(is_female=sex=='female', is_white=race=='white', is_dual=buyin_months>0) %>%
    left_join(regions, by='state') %>%
    mutate(year=year)

  # get the PDE data. record the number of claims (and days) for each bene
  # and drug
  pde = read_tsv(sprintf('../pde_%i.tsv', year)) %>%
    select(bene_id, antibiotic, days_supply) %>%
    group_by(bene_id, antibiotic) %>%
    summarize(n_claims=n(), days_supply=sum(days_supply)) %>%
    ungroup() %>%
    mutate(year=year)

  # get total numbers of PDEs
  n_claims = pde %>% group_by(bene_id) %>% summarize(n_claims=sum(n_claims))
  bene %<>%
    left_join(n_claims, by='bene_id') %>%
    replace_na(list(n_claims=0))

  # count chronic conditions. 'n_cc' means one-year cc's
  cc = read_feather(sprintf('../cc_%i.feather', year)) %>%
    rowwise() %>%
    mutate_(n_cc=sum_over(names_1yr_ccs)) %>%
    ungroup()

  bene %<>% left_join(cc, by='bene_id')

  list(bene=bene, pde=pde)
}

dat = lapply(2011:2014, load_data)
bene = lapply(dat, function(df) df$bene) %>% bind_rows
pde = lapply(dat, function(df) df$pde) %>% bind_rows
rm(dat)

# specify the constant cohort
bene %<>%
  group_by(bene_id) %>%
  mutate(n_years=n(),
         start_age=min(age),
         in_cohort=n_years==4) %>%
  ungroup()

# and include that information with the PDEs
pde %<>%
  left_join(select(bene, bene_id, in_cohort), by='bene_id')

output_table = function(x, base) write_tsv(x, sprintf('tables/tbl_%s.tsv', base))
sem = function(x) sd(x) / length(x)

summarize_totals = function(df) {
  summarize(df, n_bene=n(),
            mean_age=mean(age),
            sem_age=sem(age),
            mean_n_cc=mean(n_cc),
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
model_f(lm_cohort_bene, n_claims ~ year + start_age + is_female + is_dual + is_white + region, 'cohort_model_test')
model_f(lm_cohort_bene, n_claims ~ year + start_age*n_cc + is_female + is_dual + is_white + region, 'cohort_model2')

# models for each individual drug
subtract_min = function(x) x - min(x)
single_abx_model = function(bene, abx, frmla) {
  pde %>%
    rename(y=n_claims) %>%
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

f = function(a) single_abx_model(lm_cohort_bene, a, y ~ year + start_age*n_cc + is_female + is_dual + is_white + region)
lapply(top_abx, f) %>%
  bind_rows() %T>%
  output_table('cohort_model_abx')

# special levo & azithro models
single_abx_model(lm_cohort_bene, 'azithromycin', y ~ year + start_age + is_female + is_dual + is_white + region + COPD + ASTHMA) %>%
  output_table('cohort_model_azithro')
single_abx_model(lm_cohort_bene, 'levofloxacin', y ~ year + start_age + is_female + is_dual + is_white + region + COPD + ASTHMA) %>%
  output_table('cohort_model_levo')

#single_condition_model = function(condition_code) {
#  filter(bene, age >= 67) %>%
#    rename_(.dots=setNames(condition_code, 'has_condition')) %>%
#    glm(n_claims ~ age + is_female + has_condition + year,
#        family=quasipoisson(link='identity'), data=.) %>%
#    tidy() %>%
#    mutate(condition=condition_code)
#}

#conditions = c("AMI", "ALZHDMTA", "ATRIALFB", "CATARACT", "CHRNKIDN", "COPD", "CHF", "DIABETES", "GLAUCOMA", "HIPFRAC", "ISCHMCHT", "DEPRESSN", "OSTEOPRS", "RA_OA", "STRKETIA", "CNCRBRST", "CNCRCLRC", "CNCRPRST", "CNCRLUNG", "CNCRENDM", "ANEMIA", "ASTHMA", "HYPERL", "HYPERP", "HYPERT", "HYPOTH")
#single_condition_models = lapply(conditions, function(cnd) single_condition_model(cnd)) %>%
#  bind_rows() %T>%
#  output_table('single_condition_models')

#all_condition_formula = formula(paste0('n_claims ~ year + age + is_female + is_dual + is_white + region + ', paste0(conditions, collapse=' + ')))

# do a linear model to guess the coefficients
#all_condition_model = filter(bene, age >= 67) %>%
#  lm(all_condition_formula, data=.) %>%
#  tidy %T>%
#  output_table('model_all_condition')
