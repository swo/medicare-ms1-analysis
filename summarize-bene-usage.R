#!/usr/bin/env Rscript

# create summary beneficiary x abx file

codes = read_tsv('../db/proc/code-map.txt', col_names=c('ndc', 'abx'))

all_pde_counts = read_tsv('../data/pde_freq_2011.txt') %>%
  rename(bene=BENE_ID, all_pde=COUNT)

load.abxpde = function(fn) {
  read_tsv(fn) %>%
    rename(bene=BENE_ID, ndc=PRDSRVID, days=DAYSSPLY) %>%
    mutate(days=as.integer(days)) %>%
    left_join(codes, by='ndc') %>%
    filter(!is.na(abx)) %>%
    select(bene, abx, days)
}

county_codes = read_tsv("../db/state-county-codes/county-codes.tsv")
zipcodes = read_tsv("../db/zipcode/zipcode.tsv")
race_codes = read_tsv("../db/race-codes/race.tsv")
chronic = read_tsv("cc2011.tsv")

load.bene = function(fn, n_max=Inf) {
  # get and rearrange raw data
  read_tsv(fn, n_max=n_max) %>%
    # join state and county codes to merge with SSA list
    unite(state_county_code, STATE_CD, CNTY_CD, sep='') %>%
    # keep only the 5-digit zip
    mutate(zip=substr(BENE_ZIP, 0, 5)) %>%
    mutate(sex=as.character(factor(SEX, levels=c(0, 1, 2), labels=c('unknown', 'male', 'female')))) %>%
    select(bene=BENE_ID, age=AGE, sex, state_county_code, zip,
           race_code=RACE,
           plan_coverage_months=PLNCOVMO) %>%
    # merge (and drop) state+county code
    left_join(county_codes, by='state_county_code') %>%
    select(-state_county_code) %>%
    # merge (and drop) race codes
    left_join(race_codes, by='race_code') %>%
    select(-race_code)
}

bene = load.bene('../data/bene_match_2011.txt')
abxpde = load.abxpde('../data/abx_pde_2011.tsv')

# summarize information about each beneficiary
# (summarize over bene)
abxpde %<>%
  group_by(bene, abx) %>%
  summarize(n_claims=n(), days=sum(days))

# summary stats about beneficiaries
# (summarize again over abx for each bene)
bene_sum = abxpde %>%
  summarize(n_claims=sum(n_claims), days=sum(days))

# convert to a wide abx pde
wap = abxpde %>%
  ungroup() %>%
  select(-n_claims) %>%
  # spread to a wide format. use convert to ensure integers
  spread(abx, days, fill=0, convert=TRUE)

# merge in the summary information
wap %<>% left_join(bene_sum, by='bene')

# merge in the information from all the pdes
wap %<>% left_join(all_pde_counts, by='bene')

# merge in the information about the beneficiaries
wap %<>% left_join(bene, by='bene')

# merge in the information about the chronic conditions
wap %<>% left_join(chronic, by='bene')

# save this new file
write_tsv(wap, 'bene-abx-2011.tsv')
