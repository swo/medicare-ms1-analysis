#!/usr/bin/env Rscript

polish_year = function(y) {
  # bene
  bene_in_fn = sprintf('../data/bene_%i.tsv', y)
  bene_out_fn = sprintf('bene_%i.tsv', y)

  read_tsv(bene_in_fn, col_types=cols(zipcode='c')) %>%
    mutate(dual=BUYIN_MO>0) %>%
    select(bene_id=BENE_ID, state, sex, race, zipcode, age=AGE, dual) %>%
    write_tsv(bene_out_fn)

  # pde
  pde_in_fn = sprintf('../data/pde_%i.tsv', y)
  pde_out_fn = sprintf('pde_%i.tsv', y)

  read_tsv(pde_in_fn) %>%
    select(bene_id=BENE_ID, pde_id=PDE_ID, pde_date, antibiotic) %>%
    mutate(pde_date=dmy(pde_date)) %>%
    write_tsv(pde_out_fn)

  # chronic conditions
  cc_fn = sprintf('../data/cc_%i.tsv', y)
  cc_out_fn = sprintf('cc_%i.feather', y)
  read_tsv(cc_fn) %>%
    rename(bene_id=BENE_ID) %>%
    mutate_at(vars(-bene_id), function(x) x==3) %>%
    write_feather(cc_out_fn)
}

lapply(2011:2014, polish_year)
