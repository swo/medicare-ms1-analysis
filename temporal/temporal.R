#!/usr/bin/env Rscript

# create summary beneficiary x abx file

# - load map from NDC codes to abx names
# - load beneficiary data; keep bene_id -> state mapping
# - load antibiotic PDEs; keep bene, date, abx, days
# - merge abx and bene; keep date, abx, days, state
# - summarize information by drug, state, and year
# - summarize also by class?

library(lubridate)

codes = read_tsv('../../db/proc/code-map.txt', col_names=c('ndc', 'abx'))

load_abxpde = function(fn) {
  read_tsv(fn) %>%
    rename(bene=BENE_ID, ndc=PRDSRVID) %>%
    mutate(days=as.integer(DAYSSPLY), date=dmy(SRVC_DT), week=week(date)) %>%
    left_join(codes, by='ndc') %>%
    filter(!is.na(abx)) %>%
    mutate(abx=factor(abx)) %>%
    select(bene, abx, date, week, days)
}

load_bene = function(fn) {
  read_tsv(fn) %>%
    select(bene, state) %>%
    mutate(state=factor(state))
}

analyze = function(year) {
  abx_in_fn = sprintf('../../data/abx_pde_%s.tsv', year)
  bene_in_fn = sprintf('../bene-%s.tsv', year)
  out_fn = sprintf('temporal-abx-%s.tsv', year)
  
  abx = load_abxpde(abx_in_fn)
  bene = load_bene(bene_in_fn)
  
  dat = abx %>%
    left_join(bene, by='bene') %>%
    filter(!is.na(state), !is.na(week)) %>%
    group_by(week, state, abx) %>%
    summarize(n_rx=n(), n_days=sum(days)) %>%
    ungroup
  
  write_tsv(dat, out_fn)
}

analyze('2011')
analyze('2012')
analyze('2013')
analyze('2014')