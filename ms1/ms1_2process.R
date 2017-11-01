#!/usr/bin/env Rscript

output_table = function(x, base) write_tsv(x, sprintf('tables/%s.tsv', base))
sem = function(x) sd(x) / sqrt(length(x))

regions = read_tsv('db/census-regions.tsv') %>%
  select(state, region)

bene = read_tsv('raw_data/bene.tsv') %>%
  rename(bene_id=BENE_ID) %>%
  left_join(regions, by='state')

pde = read_tsv('raw_data/bene_pde.tsv') %>%
  rename(bene_id=BENE_ID) %>%
  rename(overall=n_pde) %>%
  left_join(bene, by=c('year', 'bene_id')) %>%
  mutate(age_group=case_when(between(.$age, 65, 75) ~ 'age65_75',
                             between(.$age, 76, 85) ~ 'age76_85',
                             between(.$age, 86, 95) ~ 'age86_95',
                             .$age > 95 ~ 'age96_'))

# summary characteristics
pop_chars = bene %>%
  group_by(year) %>%
  summarize(n_bene=n(),
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
            n_region_midwest=sum(region=='Midwest')) %T>%
  output_table('pop_chars')

n_unique_bene = length(unique(bene$bene_id)) %T>%
  write('tables/n_unique_bene.txt')

claims_by_abx = pde %>%
  select(year, overall:clindamycin) %>%
  mutate(n_bene=1) %>%
  group_by(year) %>%
  summarize_if(is.numeric, sum) %>%
  mutate_at(vars(overall:clindamycin), function(x) x/.$n_bene*1000) %>%
  select(-n_bene) %>%
  gather('abx', 'cpkp', -year) %T>%
  output_table('claims_by_abx')

# from here on, the denominators are taken from the beneficiary data grouped
# in the same way as the consumption data, so we can use the summarize_by function
summarize_by = function(df, by) {
  group_by(df, year, !!(rlang::sym(by))) %>%
    summarize(n_bene=n(), overall=sum(overall)) %>%
    ungroup() %>%
    mutate(cpkp=overall*1000/n_bene) %>%
    output_table(paste0('claims_by_', by))
}

summarize_by(pde, 'sex')
summarize_by(pde, 'region')
summarize_by(pde, 'race')
summarize_by(pde, 'age_group')

# run a model and tidy
sandwich_tidy = function(m) tidy(lmtest::coeftest(m, vcov=sandwich::vcovHC(m, type='HC3')))
model_f = function(df, frmla) {
  glm(frmla, family='poisson', data=df) %>%
    sandwich_tidy()
}

# trends by drug
drug_trend_f = function(abx) {
  frmla = str_interp("${abx} ~ year + age + n_cc + sex + race + dual + region") %>%
    as.formula
  model_f(pde, frmla) %>%
    mutate(abx=abx)
}

top_abx = c('overall', 'azithromycin', 'ciprofloxacin', 'amoxicillin',
            'cephalexin', 'tmpsmx', 'levofloxacin', 'amoxclav',
            'doxycycline', 'nitrofurantoin', 'clindamycin')

trend_models = lapply(top_abx, drug_trend_f) %>%
  bind_rows() %T>%
  output_table('trend_models')

# model overall consumption, but look at individual populations
population_models = function(term, drop_term=NULL) {
  # prepare the regression formula without that term
  covariates=c('age', 'sex', 'race', 'region')
  if (is.null(drop_term)) drop_term = term
  remaining_covariates = covariates[-which(covariates == drop_term)] %>% str_c(collapse=' + ')
  frmla = str_interp("overall ~ year + n_cc + dual + ${remaining_covariates}") %>% as.formula
  
  # get the populations
  term_values = unique(pde[[term]])
  
  lapply(term_values, function(value) {
    pde[pde[[term]]==value, ] %>%
      model_f(frmla) %>%
      mutate(population=str_c(term, '_', value))
  }) %>% bind_rows()
}

bind_rows(
  population_models('age_group', drop_term='age'),
  population_models('sex'),
  population_models('race'),
  population_models('region')
) %>%
  output_table('demography_models')
