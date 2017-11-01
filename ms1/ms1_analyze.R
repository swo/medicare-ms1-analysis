#!/usr/bin/env Rscript

output_table = function(x, base) write_tsv(x, sprintf('tables/%s.tsv', base))
sem = function(x) sd(x) / sqrt(length(x))

regions = read_tsv('db/census-regions.tsv') %>%
  select(state, region)

bene = read_tsv('raw_data/bene.tsv') %>%
  rename(bene_id=BENE_ID) %>%
  left_join(regions, by='state')

n_unique_bene = length(unique(bene$bene_id)) %T>%
  write('tables/n_unique_bene.txt')

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
glm_f = function(df, frmla, ...) eval(substitute(function(df2) glm(frmla, data=df2, ...)))(df)
sandwich_tidy = function(m) tidy(lmtest::coeftest(m, vcov=sandwich::vcovHC(m, type='HC3')))
model_f = function(df, frmla) {
  glm_f(df, frmla, family='poisson') %>%
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
population_models = function(group_type, covariate=NULL) {
  # prepare the regression formula without that term
  covariates=c('age', 'sex', 'race', 'region')
  if (is.null(covariate)) covariate = group_type
  remaining_covariates = covariates[-which(covariates == covariate)] %>% str_c(collapse=' + ')
  frmla = str_interp("overall ~ year + n_cc + dual + ${remaining_covariates}") %>% as.formula
  
  # get the populations
  pde %>%
    rename(group=!!rlang::sym(group_type)) %>%
    group_by(group) %>%
    do(model_f(., frmla)) %>%
    ungroup() %>%
    mutate(group_type=group_type)
}

bind_rows(
  population_models('age_group', covariate='age'),
  population_models('sex'),
  population_models('race'),
  population_models('region')
) %>%
  output_table('demography_models')

# appropriateness trend
pde_inapprop = read_tsv('raw_data/bene_inapprop.tsv') %>%
  filter(n_firstfill>0) %>%
  mutate(risk=n_inapprop_firstfill/n_firstfill) %>%
  rename(bene_id=BENE_ID) %>%
  left_join(bene, by=c('year', 'bene_id'))

inapprop_risk_2011 = pde_inapprop %>%
  filter(year==2011) %$%
  mean(risk) %T>%
  write('tables/inapprop_risk_2011.txt')

inapprop_trend = pde_inapprop %>%
  glm_f(risk ~ year + age + n_cc + sex + race + dual + region, weights=n_firstfill, family='binomial') %>%
  sandwich_tidy() %T>%
  output_table('inapprop_trend')

# trends in prescribing practice
dx_to_pde = read_tsv('raw_data/dx_to_pde.tsv') %>%
  rename(bene_id=BENE_ID) %>%
  left_join(bene, by=c('year', 'bene_id'))

dx_to_pde_table = dx_to_pde %>%
  select(year, dx_cat, azithromycin:clindamycin) %>%
  mutate(n_encounters=1) %>%
  group_by(year, dx_cat) %>%
  summarize_all(sum) %>%
  ungroup() %>%
  mutate_at(vars(azithromycin:clindamycin), function(x) x / .$n_encounters) %T>%
  output_table('dx_to_pde_table')

# swo: don't actually do all these combinations
dxs = c('gi_t3', 'uti_t3', 'uri', 'other_resp', 'gi',
        'asthma', 'ssti_t3', 'uti', 'ssti', 'bronchitis', 'misc_t3', 'om',
        'pharyngitis', 'pneumonia', 'misc_bacterial', 'sinusitis', 'flu',
        'om_t3', 'acne', 'viral_flu')

abx = c('azithromycin', 'ciprofloxacin', 'amoxicillin', 'cephalexin',
            'tmpsmx', 'levofloxacin', 'amoxclav',
            'doxycycline', 'nitrofurantoin', 'clindamycin')

dx_to_pde_trend_f = function (a, dx) {
  frmla = str_interp("${a} ~ year + age + n_cc + sex + race + dual + region") %>%
    as.formula
  
  now = format(Sys.time(), "%H:%M.%S")
  cat(str_interp("${now} -- abx: ${a} dx: ${dx}\n"))
  
  dx_to_pde %>%
    filter(dx_cat==dx) %>%
    glm_f(frmla, family='binomial') %>%
    sandwich_tidy()
}

dx_to_pde_trends = crossing(abx=abx, dx_cat=dxs) %>%
  group_by_all() %>%
  do(dx_to_pde_trend_f(.$abx, .$dx_cat)) %>%
  ungroup() %T>%
  output_table('dx_to_pde_trends')
