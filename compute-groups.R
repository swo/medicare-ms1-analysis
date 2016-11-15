#!/usr/bin/env Rscript

# number of data entries to read
n_max = 1e3

codes = read_tsv('../db/proc/code-map.txt', col_names=c('ndc', 'abx'))

county_codes = read_tsv("../db/state-county-codes/county-codes.tsv")
race_codes = read_tsv("../db/race-codes/race.tsv")
chronic = read_tsv("cc2011.tsv")

read_tsv('../data/abx_pde_2011.tsv') %>%
  rename(bene=BENE_ID, days=DAYSSPLY) %>%
  mutate(days=as.integer(days)) %>%
  group_by(bene) %>%
  summarize(total_days=sum(days))
  select(bene, total_days)

load.bene = function(fn, n_max=Inf) {
  # get and rearrange raw data
  read_tsv(fn, n_max=n_max) %>%
    # join state and county codes to merge with SSA list
    unite(state_county_code, STATE_CD, CNTY_CD, sep='') %>%
    # keep only the 5-digit zip
    mutate(zip=substr(BENE_ZIP, 0, 5)) %>%
    select(bene=BENE_ID, age=AGE, state_county_code, zip,
           race_code=RACE, sex,
           plan_coverage_months=PLNCOVMO) %>%
    # merge (and drop) state+county code
    left_join(county_codes, by='state_county_code') %>%
    select(-state_county_code) %>%
    # merge (and drop) race codes
    left_join(race_codes, by='race_code') %>%
    select(-race_code)
}

bene = load.bene('../data/bene_match_2011.txt', n_max=n_max)
abxpde = load.abxpde()

usage_f = function(x) ifelse(x == 0, 'none', ifelse(x <= 50, 'low', ifelse(x <=500, 'medium', 'high')))

abxpde %>%
  left_join(bene, by='bene') %>%
  left_join(chronic, by='bene') %>%
  mutate(usage=usage_f(total_days)) %>%
  group_by(usage) %>%
  summarize(n_people=n(), mean_age=mean(age), sd_age=sd(age),
            mean_com=mean(comorbidity), sd_com=sd(comorbidity),
            

# convert to a wide abx pde data frame
proto_wap = abxpde %>%
  group_by(bene, abx) %>%
  summarize(n_claims=n(), days=sum(days))

wap = proto_wap %>%
  ungroup %>%
  spread(abx, days, fill=0)

wap$total_abx = wap %>% select(-bene, -n_claims) %>% rowSums

# merge in the information from all the pdes
wap %<>% inner_join(all_pde_counts, by='bene')

# merge in the information about the beneficiaries
wap %<>% inner_join(bene, by='bene')

# merge in the information about the chronic conditions
wap %<>% inner_join(chronic, by='bene')

# save this new file
write_tsv(wap, 'bene-abx-2011.tsv')

