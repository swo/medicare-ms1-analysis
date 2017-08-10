# Consumption data
ineq = read_tsv('ineq.tsv') %>%
  group_by(drug_group, unit_type, unit) %>%
  summarize_at(vars(mean, fnz, nzgini), mean) %>%
  ungroup()

# # Aggregate usage correlates with fraction of users
# p = ggplot(all_ineq, aes(x=mean, y=fnz, color=factor(year))) +
#   geom_point() +
#   xlab('claims per beneficiary') +
#   ylab('fraction of beneficiaries w/ >= 1 claim') +
#   theme_minimal()

# # Inequality among users does not correlate well
# p = ggplot(all_ineq, aes(x=mean, y=nzgini, color=factor(year))) +
#   geom_point() +
#   xlab('claims per beneficiary') +
#   ylab('Gini coefficient among benes. w/ >= 1 claim') +
#   theme_minimal()

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

raw_abg = read_tsv('../../../antibiogram/data/abg.tsv', col_types=cols(zipcode='c')) %>%
  mutate(drug=tolower(drug), percent_nonsusceptible=100-percent_susceptible) %>%
  left_join(hrr, by=c('zipcode', 'state')) %>%
  right_join(resistance_groups, by='drug') %>%
  right_join(bug_names, by='bug') %>%
  select(bug=bug_short, drug_group, state, hrr, percent_nonsusceptible, n_isolates)

# Keep only the bug/drug combos that appear in enough states
good_combos = raw_abg %>%
  semi_join(possible_bug_drug, by=c('bug', 'drug_group')) %>%
  group_by(bug, drug_group) %>%
  summarize_at(vars(state, hrr), function(x) length(unique(x))) %>%
  ungroup() %>%
  filter(state >= 35) %>%
  select(bug, drug_group)

abg = raw_abg %>%
  semi_join(good_combos, by=c('bug', 'drug_group')) %>%
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
    group_by_(df, 'bug', 'drug_group', place_name) %>%
    summarize(mean_percent_nonsusceptible=weighted.mean(percent_nonsusceptible, sqrt(n_isolates)),
              n_place_antibiograms=n()) %>%
    mutate(place_weight=sqrt(n_place_antibiograms)) %>%
    ungroup()
}

state_abg = summarize_abg(abg, 'state')
hrr_abg = summarize_abg(abg, 'hrr')

linear_model = function(df, y, xs) {
  frmla = as.formula(str_interp("${y} ~ ${str_c(xs, collapse=' + ')}"))
  lm(formula=frmla, weights=place_weight, data=df) %>%
    tidy %>%
    mutate(ci=1.96*std.error,
           ci_low=estimate-ci,
           ci_high=estimate+ci) %>%
    select(term, estimate, ci_low, ci_high, p.value)
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

models = function(df, y, uni_x, multi_xs) {
  bind_rows(
    linear_model(df, y, uni_x) %>% mutate(model='univariate'),
    spearman_model(df, y, uni_x) %>% mutate(model='spearman'),
    linear_model(df, y, multi_xs) %>% mutate(model='multivariate')
  ) %>%
    mutate(n_data=nrow(df))
}

hrr_results = ineq %>%
  filter(unit_type=='hrr') %>%
  mutate(hrr=as.integer(unit)) %>%
  right_join(hrr_abg, by=c('hrr', 'drug_group')) %>%
  group_by(bug, drug_group) %>%
  do(models(., 'mean_percent_nonsusceptible', 'mean', c('fnz', 'nzgini'))) %>%
  ungroup()
  
state_results = ineq %>%
  filter(unit_type=='state') %>%
  rename(state=unit) %>%
  right_join(state_abg, by=c('state', 'drug_group')) %>%
  group_by(bug, drug_group) %>%
  do(models(., 'mean_percent_nonsusceptible', 'mean', c('fnz', 'nzgini'))) %>%
  ungroup()

bind_rows(
  hrr_results %>% mutate(unit_type='hrr'),
  state_results %>% mutate(unit_type='state')
) %>% write_tsv('results.tsv')
