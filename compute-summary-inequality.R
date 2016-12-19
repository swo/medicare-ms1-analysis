#!/usr/bin/env Rscript

usage_f = function(x) {
  out = rep(0, length.out=length(x))
  out[x > 0] = ntile(x[x > 0], 3)
  out
}

x = bene %>%
  filter(age >= 65) %>%
  left_join(bene_sum, by='bene') %>%
  left_join(chronic, by='bene') %>%
  mutate(days=ifelse(is.na(days), 0, days)) %>%
  mutate(is_white=ifelse(race=='white', 1, 0),
         is_male=ifelse(sex=='male', 1, 0))