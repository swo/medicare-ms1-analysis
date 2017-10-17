#!/usr/bin/env Rscript

library(rstan)

my_abx = 'quinolone'
my_bug = 'E. coli'

# Resistance data
resistance_groups = read_tsv('../resistance_groups.tsv')
bug_names = read_tsv('../bug_names.tsv')
possible_bug_drug = read_tsv('../bug_drug.tsv')

hrr = read_tsv('../../../db/hrr/hrr.tsv', col_types=list(hrr='c')) %>%
  filter(year==2014) %>%
  select(zipcode, hrr) %>%
  left_join(zipcodes, by='zipcode') %>%
  select(-state)

abg = read_tsv('../../../../antibiogram/data/abg.tsv', col_types=cols(zipcode='c')) %>%
  mutate(drug=tolower(drug), percent_nonsusceptible=100-percent_susceptible) %>%
  left_join(hrr, by='zipcode') %>%
  right_join(resistance_groups, by='drug') %>%
  right_join(bug_names, by='bug') %>%
  select(bug=bug_short, drug_group, hrr, percent_nonsusceptible, raw_n_isolates=n_isolates) %>%
  # get just the data I want
  filter(drug_group==my_abx, bug==my_bug) %>%
  # fix any NAs with medians
  mutate(n_isolates=if_else(is.na(raw_n_isolates) | raw_n_isolates==0,
                            as.integer(median(raw_n_isolates[!is.na(raw_n_isolates) & raw_n_isolates > 0])),
                            raw_n_isolates)) %>%
  mutate(n_resistant=round(percent_nonsusceptible/100*n_isolates)) %>%
  arrange(hrr)

Size = abg %>% count(hrr) %>% pull(n)
Iso = abg$n_isolates
Res = abg$n_resistant
A = nrow(abg)

states = unique(abg$hrr)
S = length(states)

Cons = ineq %>%
  filter(drug_group==my_abx, unit_type=='hrr', unit %in% states) %T>%
  { stopifnot(all(.$unit==states)) } %$%
  mean

fit = stan(file='model.stan',
           data=c('S', 'A', 'Size', 'Iso', 'Res', 'Cons'),
           iter=1000,
           chains=2)

post_vals = extract(fit, pars=c('beta0', 'beta1', 'psi')) %>%
  gather('key', 'value', everything()) %>%
  group_by(key) %>%
  summarize(mean=mean(value), sd=sd(value), n=n())