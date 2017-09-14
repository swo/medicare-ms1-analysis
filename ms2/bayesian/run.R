#!/usr/bin/env Rscript

library(rstan)

my_abx = 'tmp_smx'
my_bug = 'E. coli'

# Resistance data
resistance_groups = read_tsv('../resistance_groups.tsv')
bug_names = read_tsv('../bug_names.tsv')
possible_bug_drug = read_tsv('../bug_drug.tsv')

abg = read_tsv('../../../../antibiogram/data/abg.tsv', col_types=cols(zipcode='c')) %>%
  mutate(drug=tolower(drug), percent_nonsusceptible=100-percent_susceptible) %>%
  right_join(resistance_groups, by='drug') %>%
  right_join(bug_names, by='bug') %>%
  select(bug=bug_short, drug_group, state, percent_nonsusceptible, raw_n_isolates=n_isolates) %>%
  # get just the data I want
  filter(drug_group==my_abx, bug==my_bug) %>%
  # fix any NAs with medians
  mutate(n_isolates=if_else(is.na(raw_n_isolates) | raw_n_isolates==0,
                            as.integer(median(raw_n_isolates[!is.na(raw_n_isolates) & raw_n_isolates > 0])),
                            raw_n_isolates)) %>%
  mutate(n_resistant=round(percent_nonsusceptible/100*n_isolates)) %>%
  arrange(state)

Size = abg %>% count(state) %>% pull(n)
Iso = abg$n_isolates
Res = abg$n_resistant
A = nrow(abg)

# Consumption data
ineq = read_tsv('../ineq.tsv') %>%
  filter(drug_group==my_abx) %>%
  filter(unit_type=='state') %>%
  rename(state=unit) %>%
  filter(state %in% abg$state) %>%
  group_by(state) %>%
  mutate(f=n_bene/sum(n_bene)) %>%
  ungroup()

NC = max(ineq$n_claims) + 1
states = unique(ineq$state)
S = length(states)

Fmat = matrix(data=0, nrow=S, ncol=NC)

for (ineq_i in 1:nrow(ineq)) {
  i = match(ineq$state[ineq_i], states)
  j = ineq$n_claims[ineq_i] + 1
  f = ineq$f[ineq_i]
  
  Fmat[i, j] = f
}

fit = stan(file='model.stan',
           data=c('S', 'A', 'NC', 'Size', 'Iso', 'Res', 'Fmat'),
           iter=1000,
           chains=1)
