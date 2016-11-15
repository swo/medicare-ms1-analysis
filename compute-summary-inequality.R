#!/usr/bin/env Rscript

# Treating all usage as the same, what's the inequality in usage?

# read in accessory files
race_codes = read_tsv("../db/race-codes/race.tsv")

# read in the summary data
pde_summary = read_tsv('summary-abx-pde-2011.tsv')

# readin the chronic condition data
cc = read_tsv('cc2011.tsv')

usage_f = function(x) ifelse(x==0, 'none', ifelse(x <= 50, 'low', ifelse(x <= 500, 'medium', 'high')))
usage_f = function(x) {
  ordered_x = x[order(x)]
  # assume sorted
  out = rep(NA, length.out=length(x))
  last_zero_index = max(which(ordered_x == 0))
  tertile_size = ((length(x) - last_zero_index)) / 3 %>% floor
  out[1: last_zero_index] = 0
  out[last_zero_index: (last_zero_index + tertile_size)] = 1
  out[(last_zero_index + tertile_size): (last_zero_index + 2 * tertile_size)] = 2
  out[(last_zero_index + 2 * tertile_size): (length(x) + 1)] = 3
  out[order(x)]
}
  

# read in the beneficiary data
bene = read_tsv('../data/bene_match_2011.txt') %>%
  rename(bene=BENE_ID, race_code=RACE) %>%
  mutate(sex=ifelse(SEX==1, "male", ifelse(SEX==2, "female", "unknown"))) %>%
  mutate(zip=substr(BENE_ZIP, 0, 5)) %>%
  select(bene, race_code, age=AGE, sex) %>%
  # merge (and drop) race codes
  left_join(race_codes, by='race_code') %>%
  select(-race_code) %>%
  # merge chronic conditions
  left_join(cc, by='bene') %>%
  left_join(pde_summary, by='bene') %>%
  mutate_at(vars(matches('total_')), function(x) ifelse(is.na(x), 0, x)) %>%
  mutate(usage=usage_f(total_days))

write_tsv(bene, 'tmp')