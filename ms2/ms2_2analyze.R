# Consumption data
ineq = read_tsv('ineq.tsv')

# regions = read_tsv('../../db/census-regions/census-regions.tsv') %>%
#   select(state, state_abbreviation, region)
# 
# zipcodes = read_tsv('../../db/zipcode/zipcode.tsv') %>%
#   select(zipcode=zip, state)
# 
# read_years = function(base_fn, years) {
#   lapply(years, function(y) {
#       read_tsv(sprintf(base_fn, y)) %>% mutate(year=y)
#     }) %>%
#     bind_rows
# }
# 
# hrr = read_tsv('../../db/hrr/hrr.tsv') %>%
#   select(zipcode, hrr) %>%
#   left_join(zipcodes, by='zipcode') %>%
#   left_join(regions, by='state')
# 
# years = 2011:2014
# all_ineq = read_years('data/state_ineq_all_%i.tsv', years) %>%
#   left_join(regions, by='state')
# state_ineq = read_years('data/state_ineq_%i.tsv', years)
# hrr_ineq = read_years('data/hrr_ineq_%i.tsv', years)

# Aggregate usage correlates with fraction of users
p = ggplot(all_ineq, aes(x=mean, y=fnz, color=factor(year))) +
  geom_point() +
  xlab('claims per beneficiary') +
  ylab('fraction of beneficiaries w/ >= 1 claim') +
  theme_minimal()

ggsave('table/fig_all_ineq.pdf', plot=p)

# Inequality among users does not correlate well
p = ggplot(all_ineq, aes(x=mean, y=nzgini, color=factor(year))) +
  geom_point() +
  xlab('claims per beneficiary') +
  ylab('Gini coefficient among benes. w/ >= 1 claim') +
  theme_minimal()

ggsave('table/fig_all_gini.pdf', plot=p)

# Resistance data

resistance_groups = read_tsv('resistance_groups.tsv')
bug_names = read_tsv('bug_names.tsv')

abg = read_tsv('../../../antibiogram/data/abg.tsv', col_types=cols(zipcode='c')) %>%
  mutate(drug=tolower(drug), percent_nonsusceptible=100-percent_susceptible) %>%
  left_join(hrr, by=c('zipcode', 'state')) %>%
  right_join(resistance_groups, by='drug') %>%
  select(bug, drug_group, state, hrr, percent_nonsusceptible, n_isolates)

possible_bug_drug = read_tsv('bug_drug.tsv')

# Keep only bug/drug combinations that are in at least 75 HRRs
min_n_states = 35
good_combos = semi_join(abg, possible_bug_drug, by=c('bug', 'drug_group')) %>%
  select(bug, drug_group, state) %>%
  distinct %>%
  count(bug, drug_group) %>% rename(n_places=n) %>%
  filter(n_places >= min_n_states)

abg %<>% right_join(good_combos, by=c('bug', 'drug_group')) %>%
  rename(bug_long=bug) %>% left_join(bug_names, by='bug_long') %>% select(-bug_long) %>%
  filter(!is.na(percent_nonsusceptible), !is.na(n_isolates))

# check that all the short names made it in
if (any(is.na(abg$bug))) stop('missing bug names')

summarize_abg = function(df, place_name) {
    group_by_(df, 'bug', 'drug_group', place_name) %>%
    summarize(mean_percent_nonsusceptible=weighted.mean(percent_nonsusceptible, sqrt(n_isolates)),
              n_place_antibiograms=n()) %>%
    mutate(place_weight=sqrt(n_place_antibiograms))
}

state_abg = summarize_abg(abg, 'state')
hrr_abg = summarize_abg(abg, 'hrr')

linear_model = function(df, dependent_var, explanatory_vars) {
  f = paste0(dependent_var, ' ~ ', paste0(explanatory_vars, collapse = ' + ')) %>% formula
  lm(formula=f, weights=place_weight, data=df) %>%
    tidy %>%
    select(term, estimate, p.value)
}

spearman_model = function(df, dependent_var, explanatory_var) {
  cor.test(df[[dependent_var]], df[[explanatory_var]], method='spearman') %>%
    tidy %>%
    mutate(term=explanatory_var, method='spearman') %>%
    select(term, estimate, p.value)
}

multivariate_term = 'nzgini'
models = function(x) {
  dep = 'mean_percent_nonsusceptible'
  bind_rows(
    linear_model(x, dep, 'mean') %>% mutate(model='univariate'),
    linear_model(x, dep, c('fnz', multivariate_term)) %>% mutate(model='multivariate'),
    spearman_model(x, dep, 'mean') %>% mutate(model='spearman')
  ) %>%
    mutate(n_data=nrow(x))
}

results = left_join(state_abg, state_ineq, by=c('state', 'drug_group')) %>%
  group_by(year, bug, drug_group) %>%
  do(models(.)) %>%
  ungroup

state_ineq %>%
  group_by(state, drug_group) %>%
  summarize_at(vars(mean, fnz, nzgini), mean) %>%
  ungroup() %>%
  right_join(state_abg, by=c('state', 'drug_group')) %>%
  group_by(bug, drug_group) %>%
  do(models(.)) %>%
  ungroup %>%
  write_tsv('table/state_average_result.tsv')

hrr_results = left_join(hrr_abg, hrr_ineq, by=c('hrr', 'drug_group')) %>%
  group_by(year, bug, drug_group) %>%
  do(models(.)) %>%
  ungroup

write_tsv(results, 'table/state_result.tsv')
write_tsv(hrr_results, 'table/hrr_result.tsv')
