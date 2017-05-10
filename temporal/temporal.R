#!/usr/bin/env Rscript

library(lubridate)

summarize_weeks = function(x) {
  mutate(x, week=week(service_date)) %>%
    group_by(state, week, antibiotic) %>%
    summarize(n_claims=n(), n_days=sum(days_supply)) %>%
    ungroup
}

for (y in 2011:2014) {
  pde = sprintf('../pde_%i.tsv', y) %>% read_tsv
  bene = sprintf('../bene_%i.tsv', y) %>%
    read_tsv %>%
    filter(between(age, 66, 96), sex %in% c('male', 'female'), hmo_months==0)

  pde %>%
    left_join(select(bene, bene_id, state), by='bene_id') %>%
    filter(!is.na(state)) %>%
    summarize_weeks %>%
    write_tsv(sprintf('temporal_abx_%i.tsv', y))
}
