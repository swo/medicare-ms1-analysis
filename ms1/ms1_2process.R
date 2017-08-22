library(forcats)

# load data
# NB: I put DC into the South
subtract_min = function(x) x - min(x)

bene = read_tsv('data/bene.tsv') %>%
  mutate(region=factor(region, levels=c('Northeast', 'West', 'Midwest', 'South')),
         sex=factor(sex, levels=c('male', 'female')),
         race=factor(race, levels=c('white', 'black', 'Hispanic', 'other')),
         dual=buyin_months>0) %>%
  mutate(age=age-1) %>%
  mutate(year0=subtract_min(year),
         age0=subtract_min(age))

pde = read_tsv('data/pde.tsv')

dx = read_tsv('data/dx.tsv') %>%
  left_join(bene, by=c('year', 'bene_id'))

dx_pde = read_tsv('data/dx_pde.tsv')

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
            n_female=sum(sex=='female'),
            n_white=sum(race=='white'),
            n_dual=sum(dual),
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

n_unique_bene = bene$bene_id %>%
  unique %>% length %T>%
  write('tables/n_unique_bene.tsv')

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

top_abx = c('azithromycin', 'ciprofloxacin', 'amoxicillin', 'cephalexin',
            'trimethoprim/sulfamethoxazole', 'levofloxacin', 'amoxicillin/clavulanate',
            'doxycycline', 'nitrofurantoin', 'clindamycin')

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

abx_frmla = y ~ year0 + age0 + n_cc + sex + race + dual + region

# get models and predictions for all drugs combined
overall_model = bene %>%
  rename(y=n_claims) %>%
  model_f(abx_frmla, 'overall')

# model overall consumption, but look at individual populations
covariates = c('age0', 'sex', 'race', 'dual', 'region')
reduced_model = function(df, missing_term, name) {
  remaining_covariates = covariates[-which(covariates == missing_term)] %>% str_c(collapse=' + ')
  frmla = str_interp("y ~ year0 + n_cc + ${remaining_covariates}") %>% as.formula
  
  df %>%
    rename(y=n_claims) %>%
    model_f(frmla, name)
}

bind_rows(
  bene %>% filter(between(age, 65, 75)) %>% reduced_model('age0', 'age65_75'),
  bene %>% filter(between(age, 76, 85)) %>% reduced_model('age0', 'age76_85'),
  bene %>% filter(between(age, 86, 95)) %>% reduced_model('age0', 'age86_95'),
  bene %>% filter(sex=='female') %>% reduced_model('sex', 'female'),
  bene %>% filter(sex=='male') %>% reduced_model('sex', 'male'),
  bene %>% filter(race=='white') %>% reduced_model('race', 'white'),
  bene %>% filter(race=='black') %>% reduced_model('race', 'black'),
  bene %>% filter(race=='Hispanic') %>% reduced_model('race', 'Hispanic'),
  bene %>% filter(race=='other') %>% reduced_model('race', 'other'),
  bene %>% filter(region=='South') %>% reduced_model('region', 'South'),
  bene %>% filter(region=='Midwest') %>% reduced_model('region', 'Midwest'),
  bene %>% filter(region=='West') %>% reduced_model('region', 'West'),
  bene %>% filter(region=='Northeast') %>% reduced_model('region', 'Northeast')
) %>%
  output_table('model_overall_demography')

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

## DIAGNOSES ##
# total counts for each diagnosis
# (to be adjusted to diagnoses per 1k ppl)
dx_counts = dx %>%
  count(year, diagnosis_category) %>%
  ungroup() %>%
  left_join(select(totals, year, n_bene), by=c('year')) %>%
  mutate(dxpkp=1000*n/n_bene) %T>%
  output_table('dx_counts')

# adjusted rates
# unit of analysis is the individual
# Poisson regression # dx ~ year + covariates
dx_model_data = dx %>%
  count(year, bene_id, diagnosis_category) %>%
  ungroup() %>%
  rename(y=n)

dx_models = lapply(unique(dx$diagnosis_category), function(x) {
  dx_model_data %>%
    filter(diagnosis_category==x) %>%
    right_join(bene, by=c('year', 'bene_id')) %>%
    replace_na(list(y=0)) %>%
    glm(abx_frmla, data=., family='poisson') %>%
    tidy %>%
    mutate(diagnosis_category=x)
}) %>%
  bind_rows() %T>%
  output_table('dx_models')

# count delays
dx_delays = dx_pde %>%
  count(delay) %T>%
  output_table('dx_delay')

# for each drug, how often is a diagnosis present?
pde_dx_data = pde %>%
  left_join(select(dx_pde, year, pde_id, dx_id, diagnosis_category), by=c('year', 'pde_id')) %>%
  replace_na(list(diagnosis_category='none')) %>%
  left_join(bene, by=c('year', 'bene_id'))

# (NB: PDE can get double-counted, so need to use raw PDE counts
# as denominator)
pde_dx_counts = pde_dx_data %>%
  count(year, antibiotic, diagnosis_category) %>%
  ungroup() %>%
  rename(n_pde_with_dx=n) %>%
  left_join(select(claims_by_abx, year, antibiotic, n_pde_denom=n), by=c('year', 'antibiotic')) %>%
  mutate(f_pde_with_dx=n_pde_with_dx/n_pde_denom) %T>%
  output_table('pde_dx_counts')

# for each drug, how does the frequecy of a diagnosis
# change with time, when adjusted for covariates?
# SWO: don't need this?
# pde_dx_adjust_f = function(abx, dx) {
#   pde_dx_data %>%
#     filter(antibiotic==abx) %>%
#     mutate(y=diagnosis_category==dx) %>%
#     glm(abx_frmla, data=., family='binomial') %>%
#     tidy %>%
#     mutate(antibiotic=abx, diagnosis_category=dx)
# }
# 
# pde_dx_models = crossing(antibiotic=top_abx, diagnosis_category=c('none', unique(dx$diagnosis_category))) %>%
#   rowwise() %>%
#   do(pde_dx_adjust_f(.$antibiotic, .$diagnosis_category)) %T>%
#   output_table('pde_dx_model')

# for each diagnosis, how often is each drug present?
dx_pde_data = dx %>%
  left_join(select(dx_pde, year, pde_id, dx_id, antibiotic), by=c('year', 'dx_id')) %>%
  replace_na(list(antibiotic='none'))

dx_pde_counts = dx_pde_data %>%
  count(year, diagnosis_category, antibiotic) %>%
  ungroup() %>%
  rename(n_dx_with_pde=n) %>%
  left_join(select(dx_counts, year, diagnosis_category, n_dx_denom=n), by=c('year', 'diagnosis_category')) %>%
  mutate(f_dx_with_pde=n_dx_with_pde/n_dx_denom) %T>%
  output_table('dx_pde_counts')

# for each diagnosis, how does the frequency of giving a drug
# change over time, adjusted for covariates?
dx_pde_adjust_f = function(dx, abx) {
  dx_pde_data %>%
    filter(diagnosis_category==dx) %>%
    mutate(y=antibiotic==abx) %>%
    glm(abx_frmla, data=., family='binomial') %>%
    tidy %>%
    mutate(antibiotic=abx, diagnosis_category=dx)
}

top_dx = dx_counts %>%
  group_by(diagnosis_category) %>%
  summarize(dxpkp=max(dxpkp)) %>%
  arrange(desc(dxpkp)) %>%
  head(17) %$%
  diagnosis_category

dx_pde_models = crossing(diagnosis_category=top_dx, antibiotic=top_abx) %>%
  rowwise() %>%
  do(dx_pde_adjust_f(.$diagnosis_category, .$antibiotic)) %>%
  ungroup() %T>%
  output_table('dx_pde_model')

# risk ratio
#rr = function(or, p) or / (1 - p + (p * or))
