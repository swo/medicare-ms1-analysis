#!/usr/bin/env Rscript

# create summary beneficiary x abx file

# - load map from NDC codes to abx names
# - load beneficiary data; keep bene_id -> state mapping
# - load antibiotic PDEs; keep bene, date, abx, days
# - merge abx and bene; keep date, abx, days, state
# - summarize information by drug, state, and year
# - summarize also by class?

library(lubridate)

summarize_weeks = function(x) {
  mutate(x, week=week(service_date)) %>%
    group_by(state, week, antibiotic) %>%
    summarize(n_claims=n(), n_days=sum(days_supply)) %>%
    ungroup
}

for (y in 2011:2014) {
  pde = sprintf('../pde_%i.tsv', y) %>% read_tsv
  bene = sprintf('../bene_%i.tsv', y) %>% read_tsv

  pde %>%
    left_join(select(bene, bene_id, state), by='bene_id') %>%
    summarize_weeks %>%
    write_tsv(sprintf('temporal_abx_%i.tsv', y))
}