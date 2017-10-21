# load data
# NB: I put DC into the South

# bene
bene = read_tsv('data/bene.tsv') %>%
  mutate(region=factor(region, levels=c('Northeast', 'West', 'Midwest', 'South')),
         sex=factor(sex, levels=c('male', 'female')),
         race=factor(race, levels=c('white', 'black', 'Hispanic', 'other'))) %>%
  mutate(age=age-1)

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

# pde
pde = read_tsv('data/pde.tsv') %>%
  left_join(bene, by=c('year', 'bene_id'))

pde_firstfill = pde %>% filter(fill_num==0)

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
sandwich_tidy = function(m) tidy(lmtest::coeftest(m, vcov=sandwich::vcovHC(m, type='HC3')))
model_f = function(df, frmla, name) {
  glm_f(df, frmla, family='poisson') %>%
    sandwich_tidy() %>%
    mutate(name=name)
}

abx_frmla = y ~ year + age + n_cc + sex + race + dual + region

# get models and predictions for all drugs combined
overall_model = bene %>%
  rename(y=n_claims) %>%
  model_f(abx_frmla, 'overall')

overall_model_firstfill = bene %>%
  rename(y=n_claims_firstfill) %>%
  model_f(abx_frmla, 'overall')

# model overall consumption, but look at individual populations
covariates = c('age', 'sex', 'race', 'dual', 'region')
reduced_model = function(df, missing_term, name) {
  remaining_covariates = covariates[-which(covariates == missing_term)] %>% str_c(collapse=' + ')
  frmla = str_interp("y ~ year + n_cc + ${remaining_covariates}") %>% as.formula

  df %>%
    rename(y=n_claims) %>%
    model_f(frmla, name)
}

bind_rows(
  bene %>% filter(between(age, 65, 75)) %>% reduced_model('age', 'age65_75'),
  bene %>% filter(between(age, 76, 85)) %>% reduced_model('age', 'age76_85'),
  bene %>% filter(between(age, 86, 95)) %>% reduced_model('age', 'age86_95'),
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
single_f = function(pde, abx, frmla) {
  pde %>%
    filter(antibiotic==abx) %>%
    count(year, bene_id) %>% ungroup() %>% rename(y=n) %>%
    right_join(bene, by=c('year', 'bene_id')) %>%
    replace_na(list(y=0)) %>%
    model_f(frmla, abx)
}

single_models = lapply(top_abx, function(a) single_f(pde, a, abx_frmla))
single_models_firstfill = lapply(top_abx, function(a) single_f(pde_firstfill, a, abx_frmla))

# report the model parameters
bind_rows(overall_model, single_models) %>%
  output_table('model_abx')

bind_rows(overall_model_firstfill, single_models_firstfill) %>%
  output_table('model_abx_firstfill')
