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
  group_by(year, unit_type, unit, drug_group) %>%
  summarize(mean=weighted.mean(n_claims, w=n_bene),
            fnz=sum(n_bene[n_claims > 0]) / sum(n_bene),
            mup=mean/fnz) %>%
  group_by(unit_type, unit, drug_group) %>%
  summarize_at(vars(mean, fnz, mup), mean) %>%
  ungroup() %>%
  # get just the data I want
  filter(unit_type=='state', drug_group==my_abx) %>%
  rename(state=unit) %>%
  filter(state %in% abg$state)

Cons = ineq$mean
Fnz = ineq$fnz
Mup = ineq$mup
S = nrow(ineq)

# check that we got the Sizes all correct
pos = 0;
for (i in 1:S) {
  for (j in 1:Size[i]) {
    stopifnot(abg$state[pos + j] == ineq$state[i]);
  }
  pos = pos + Size[i];
}

fit = stan(file='probs3.stan',
           data=c('S', 'A', 'Size', 'Iso', 'Res', 'Fnz', 'Mup'),
           iter=1000,
           chains=1)

lm_fit = abg %>%
  group_by(state) %>%
  summarize(y=weighted.mean(percent_nonsusceptible/100, w=sqrt(n_isolates)),
            n_abg=n()) %>%
  left_join(ineq, by='state') %>%
  lm(y ~ fnz + mup, data=., weights=n_abg)
