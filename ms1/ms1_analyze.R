#!/usr/bin/env Rscript

output_table = function(x, base) write_tsv(x, sprintf('tables/%s.tsv', base))

bene = read_tsv('raw_data/bene.tsv') %>%
  rename(overall=n_pde,
         overall_inapprop=n_pde_inapprop,
         overall_wrefill=n_pde_wrefill) %>%
  mutate(overall_approp=overall - overall_inapprop) %>%
  mutate(age_group=case_when(between(.$age, 65, 75) ~ 'age65_75',
                             between(.$age, 76, 85) ~ 'age76_85',
                             between(.$age, 86, 95) ~ 'age86_95',
                             .$age > 95 ~ 'age96_'))

n_unique_bene = length(unique(bene$BENE_ID)) %T>%
  write('tables/n_unique_bene.txt')

# summary characteristics
pop_chars = bene %>%
  group_by(year) %>%
  summarize(n_bene=n(),
            mean_age=mean(age),
            sd_age=sd(age),
            mean_n_cc=mean(n_cc),
            sd_n_cc=sd(n_cc),
            n_female=sum(sex=='female'),
            n_white=sum(race=='white'),
            n_dual=sum(dual),
            n_region_south=sum(region=='South'),
            n_region_west=sum(region=='West'),
            n_region_northeast=sum(region=='Northeast'),
            n_region_midwest=sum(region=='Midwest')) %T>%
  output_table('pop_chars')

claims_by_fill = bene %>%
  select(overall, overall_wrefill) %>%
  summarize_all(sum) %>%
  mutate(refills=overall_wrefill - overall) %T>%
  output_table('claims_by_fill')

claims_by_abx = bene %>%
  select(year, overall, overall_inapprop, overall_approp, azithromycin:clindamycin) %>%
  mutate(n_bene=1) %>%
  group_by(year) %>%
  summarize_if(is.numeric, sum) %>%
  mutate_at(vars(-year, -n_bene), function(x) x/.$n_bene*1000) %>%
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

summarize_by(bene, 'sex')
summarize_by(bene, 'region')
summarize_by(bene, 'race')
summarize_by(bene, 'age_group')

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
  model_f(bene, frmla) %>%
    mutate(abx=abx)
}

top_abx = c('overall', 'overall_approp', 'overall_inapprop', 'overall_wrefill',
            'azithromycin', 'ciprofloxacin', 'amoxicillin',
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
  bene %>%
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

# trends in prescribing practice
dx_to_pde = read_tsv('raw_data/dx_to_pde.tsv')

dx_to_pde_table = dx_to_pde %>%
  select(year, dx_cat, azithromycin:clindamycin) %>%
  mutate(n_encounters=1) %>%
  group_by(year, dx_cat) %>%
  summarize_all(sum) %>%
  ungroup() %>%
  mutate_at(vars(azithromycin:clindamycin), function(x) x / .$n_encounters) %T>%
  output_table('dx_to_pde_table')

dx_to_pde_trend_f = function (a, dx) {
  frmla = str_interp("${a} ~ year + age + n_cc + sex + race + dual + region") %>%
    as.formula
  
  # now = format(Sys.time(), "%H:%M.%S")
  # cat(str_interp("${now} -- abx: ${a} dx: ${dx}\n"))
  
  dx_to_pde %>%
    filter(dx_cat==dx) %>%
    glm_f(frmla, family='binomial') %>%
    sandwich_tidy()
}

dx_to_pde_trend_pairs = bind_rows(
  crossing(abx=c('azithromycin', 'levofloxacin', 'amoxclav'),
           dx_cat=c('pneumonia', 'sinusitis', 'uri', 'bronchitis', 'pharyngitis', 'asthma', 'other_resp')),
  crossing(abx=c('ciprofloxacin', 'tmpsmx', 'nitrofurantoin'),
           dx_cat=c('uti', 'uti_t3')),
  crossing(abx=c('cephalexin', 'tmpsmx', 'clindamycin'),
           dx_cat=c('ssti', 'ssti_t3'))
)

dx_to_pde_trends = dx_to_pde_trend_pairs %>%
  group_by_all() %>%
  do(dx_to_pde_trend_f(.$abx, .$dx_cat)) %>%
  ungroup() %T>%
  output_table('dx_to_pde_trends')
