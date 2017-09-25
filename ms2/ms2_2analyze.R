#!/usr/bin/env Rscript

# Consumption data
ineq = read_tsv('ineq.tsv') %>%
  group_by(year, unit_type, unit, drug_group) %>%
  summarize(mean=weighted.mean(n_claims, w=n_bene),
            fnz=sum(n_bene[n_claims > 0]) / sum(n_bene),
            mup=mean/fnz) %>%
  group_by(unit_type, unit, drug_group) %>%
  summarize_at(vars(mean, fnz,mup), mean) %>%
  ungroup()

# Resistance data
resistance_groups = read_tsv('resistance_groups.tsv')
bug_names = read_tsv('bug_names.tsv')
possible_bug_drug = read_tsv('bug_drug.tsv')

regions = read_tsv('../../db/census-regions/census-regions.tsv') %>%
  select(state, state_abbreviation, region)

zipcodes = read_tsv('../../db/zipcode/zipcode.tsv') %>%
  select(zipcode=zip, state)

# Just use 2014 HRRs
hrr = read_tsv('../../db/hrr/hrr.tsv') %>%
  filter(year==2014) %>%
  select(zipcode, hrr) %>%
  left_join(zipcodes, by='zipcode') %>%
  left_join(regions, by='state')

# Associate each HRR with a single region by just taking the first one
hrr_region = hrr %>%
  group_by(hrr) %>%
  summarize(region=first(region))

abg = read_tsv('../../../antibiogram/data/abg.tsv', col_types=cols(zipcode='c')) %>%
  mutate(drug=tolower(drug), percent_nonsusceptible=100-percent_susceptible) %>%
  left_join(hrr, by=c('zipcode', 'state')) %>%
  right_join(resistance_groups, by='drug') %>%
  right_join(bug_names, by='bug') %>%
  select(bug=bug_short, drug_group, state, hrr, percent_nonsusceptible, n_isolates, has_inpatient) %>%
  # filter for only the "good" combinations
  semi_join(possible_bug_drug, by=c('bug', 'drug_group')) %>%
  # compute the median number of isolates for each bug/drug combo
  (function(df) {
    filter(df, !is.na(n_isolates)) %>%
      group_by(bug, drug_group) %>%
      summarize(median_n_isolates=as.integer(median(n_isolates))) %>%
      ungroup() %>%
      right_join(df, by=c('bug', 'drug_group')) }) %>%
  rename(raw_n_isolates=n_isolates) %>%
  mutate(n_isolates=if_else(is.na(raw_n_isolates) | raw_n_isolates==0,
                            median_n_isolates, raw_n_isolates))

# check that all the short names made it in
if (any(is.na(abg$bug))) stop('missing bug names')

summarize_abg = function(df, place_name) {
  place_name_var = enquo(place_name)
  place_name_str = quo_name(place_name_var)
  
  frmla = str_interp("percent_nonsusceptible ~ ${place_name_str} + has_inpatient + 0") %>% as.formula
  
  models = df %>%
    group_by(bug, drug_group) %>%
    do(tidy(lm(frmla, data=., weights=median_n_isolates))) %>%
    mutate(!!place_name_str := str_replace(term, str_c('^', place_name_str), '')) %>%
    select(-term)
  
  abg_counts = df %>%
    group_by(bug, drug_group, !!place_name_var) %>%
    summarize(n_place_antibiograms=n())
  
  left_join(models, abg_counts, by=c('bug', 'drug_group', place_name_str))
}

state_abg = summarize_abg(abg, state) %T>% write_tsv('abg_state.tsv')
hrr_abg = abg %>% mutate(hrr=as.character(hrr)) %>% summarize_abg(hrr) %T>% write_tsv('abg_hrr.tsv')

linear_model = function(df, y, xs) {
  frmla = as.formula(str_interp("${y} ~ ${str_c(xs, collapse=' + ')}"))
  m = lm(formula=frmla, weights=n_place_antibiograms, data=df)

  anova_res = anova(m) %>%
    tidy %>%
    filter(term != 'Residuals') %>%
    select(term, anova.p.value=p.value)

  coef_res = tidy(m) %>%
    mutate(ci=1.96*std.error,
           ci_low=estimate-ci,
           ci_high=estimate+ci) %>%
    select(term, estimate, ci_low, ci_high, p.value)

  left_join(coef_res, anova_res, by='term')
}

spearman_model = function(df, y_name, x_name) {
  y = df[[y_name]]
  x = df[[x_name]]

  m = cor.test(y, x, method='spearman')

  if (length(y) > 3) {
    cis = mada::CIrho(m$estimate, length(y), level=0.95)
  } else {
    cis = c(0, 0, 1)
  }

  data_frame(term=x_name,
             estimate=m$estimate,
             ci_low=cis[2],
             ci_high=cis[3],
             p.value=m$p.value)
}

beta_model = function(df_raw, y, xs, link='logit', eps=1e-6) {
  # replace 0 and 1 with slightly different values
  df = df_raw %>%
    mutate(y=case_when(.$y==0 ~ eps,
                       .$y==1 ~ 1 - eps,
                       TRUE ~ .$y))

  frmla = as.formula(str_interp("${y} ~ ${str_c(xs, collapse=' + ')}"))
  m = betareg::betareg(formula=frmla, weights=n_place_antibiograms, data=df, link=link)

  coef_res = tidy(m) %>%
    mutate(ci=1.96*std.error,
           ci_low=estimate-ci,
           ci_high=estimate+ci) %>%
    select(term, estimate, ci_low, ci_high, p.value)

  coef_res
}

logistic_model = function(df, y, xs) {
  frmla = as.formula(str_interp("${y} ~ ${str_c(xs, collapse='+')}"))
  m = glm(formula=frmla, weights=n_place_antibiograms, data=df, family=quasibinomial)

  tidy(m) %>%
    mutate(ci=1.96*std.error,
           ci_low=estimate-ci,
           ci_high=estimate+ci) %>%
    select(term, estimate, ci_low, ci_high, p.value)
}

models = function(df) {
  bind_rows(
    spearman_model(df, 'y', 'mean') %>% mutate(model='spearman'),
    linear_model(df, 'y', 'mean') %>% mutate(model='univariate_mean'),
    linear_model(df, 'y', 'fnz') %>% mutate(model='univariate_fnz'),
    linear_model(df, 'y', 'mup') %>% mutate(model='univariate_mup'),
    linear_model(df, 'y', c('fnz', 'mup')) %>% mutate(model='multivariate_fnz_mup'),
    linear_model(df, 'y', c('mup', 'fnz')) %>% mutate(model='multivariate_mup_fnz')
  ) %>%
    mutate(n_data=nrow(df))
}

hrr_results = ineq %>%
  filter(unit_type=='hrr') %>%
  mutate(hrr=as.integer(unit)) %>%
  right_join(hrr_abg, by=c('hrr', 'drug_group')) %>%
  left_join(hrr_region, by='hrr') %>%
  mutate(y=mean_percent_nonsusceptible/100) %>%
  group_by(bug, drug_group) %>%
  do(models(.)) %>%
  ungroup()

state_results = ineq %>%
  filter(unit_type=='state') %>%
  rename(state=unit) %>%
  left_join(regions, by='state') %>%
  right_join(state_abg, by=c('state', 'drug_group')) %>%
  mutate(y=mean_percent_nonsusceptible/100) %>%
  group_by(bug, drug_group) %>%
  do(models(.)) %>%
  ungroup()

results = bind_rows(
  hrr_results %>% mutate(unit_type='hrr'),
  state_results %>% mutate(unit_type='state')
) %T>%
  write_tsv('results.tsv')
