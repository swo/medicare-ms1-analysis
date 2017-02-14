#!/usr/bin/env Rscript

library(lubridate)

# create summary beneficiary x abx file

codes = read_tsv('../db/proc/code-map.txt', col_names=c('ndc', 'abx'))

load_abxpde = function(fn, expected_year) {
  read_tsv(fn) %>%
    rename(bene=BENE_ID, service_date=SRVC_DT, ndc=PRDSRVID, days=DAYSSPLY) %>%
    mutate(service_date=dmy(service_date), days=as.integer(days)) %>%
    filter(year(service_date) == expected_year) %>%
    left_join(codes, by='ndc') %>%
    filter(!is.na(abx)) %>%
    select(bene, abx, days)
}

load_bene = function(fn, n_max=Inf) {
  # get and rearrange raw data
  read_tsv(fn, n_max=n_max) %>%
    # keep only the 5-digit zip
    mutate(zip=substr(BENE_ZIP, 0, 5)) %>%
    mutate(sex=as.character(factor(SEX, levels=c(0, 1, 2), labels=c('unknown', 'male', 'female')))) %>%
    select(bene=BENE_ID, age=AGE, sex, zip,
           race_code=RACE,
           plan_coverage_months=PLNCOVMO) %>%
    # merge by zipcode to get state
    left_join(zipcodes, by='zip') %>%
    # merge (and drop) race codes
    left_join(race_codes, by='race_code') %>%
    select(-race_code)
}

# load the database information
county_codes = read_tsv("../db/state-county-codes/county-codes.tsv")
zipcodes = read_tsv("../db/zipcode/zipcode.tsv") %>%
  select(zip, state, state_abbrev)
race_codes = read_tsv("../db/race-codes/race.tsv")

# and get the common chronic condition data
all_chronic = read_tsv("../data/cc_all.tsv", col_types='cdddd')

summarize_year = function(year) {
  bene_in_fn = sprintf('../data/bene_match_%s.txt', year)
  abx_in_fn = sprintf('../data/abx_pde_%s.tsv', year)

  bene_out_fn = sprintf('bene-%s.tsv', year)
  abx_out_fn = sprintf('abx-pde-%s.tsv', year)
  wide_out_fn = sprintf('bene-abx-%s.tsv', year)

  bene = load_bene(bene_in_fn)
  abxpde = load_abxpde(abx_in_fn, expected_year=as.integer(year))

  year_col = sprintf('y%s', year)
  chronic = all_chronic %>%
    select_('bene', year_col) %>%
    rename_(.dots=setNames(year_col, 'chronic'))

  bene %>%
    left_join(chronic, by='bene') %>%
    write_tsv(bene_out_fn)

  # summarize information about each beneficiary
  # (summarize over bene)
  abxpde %<>%
    group_by(bene, abx) %>%
    summarize(n_claims=n(), days=sum(days))

  abxpde %>% write_tsv(abx_out_fn)

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

  # merge in the information about the beneficiaries
  wap %<>% left_join(bene, by='bene')

  # merge in the information about the chronic conditions
  wap %<>% left_join(chronic, by='bene')

  # save this new file
  write_tsv(wap, wide_out_fn)
}

summarize_year('2011')
summarize_year('2012')
summarize_year('2013')
summarize_year('2014')
