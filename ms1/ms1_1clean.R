#!/usr/bin/env Rscript

# Process the bene and pde
region = read_tsv('db/census-regions.tsv') %>%
  select(state, region)

pde = read_tsv('raw_data/pde.tsv') %>%
  mutate(first_fill=fill_num==0) %>%
  select(year, bene_id, antibiotic, first_fill)

pde_denom = pde %>%
  group_by(year, bene_id) %>%
  summarize(n_claims=n(),
            n_claims_firstfill=sum(first_fill)) %>%
  ungroup()

bene = read_tsv('raw_data/bene.tsv') %>%
  left_join(region, by='state') %>%
  mutate(dual=buyin_mo > 0) %>% select(-buyin_mo) %>%
  mutate(region=factor(region, levels=c('Northeast', 'West', 'Midwest', 'South')),
         sex=factor(sex, levels=c('male', 'female')),
         race=factor(race, levels=c('white', 'black', 'Hispanic', 'other'))) %>%
  mutate(age=age-1) %>%
  left_join(pde_denom, by=c('year', 'bene_id')) %>%
  replace_na(list(n_claims=0, n_claims_firstfill=0))

saveRDS(bene, 'data/bene.rds')
saveRDS(pde, 'data/pde.rds')
