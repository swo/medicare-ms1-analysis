#!/usr/bin/env Rscript

abx_class = read_tsv('abx_class.tsv')

polish_bene = function(df) {
  rename(df, bene_id=BENE_ID, age=AGE, buyin_months=BUYIN_MO)
}

polish_pde = function(df) {
  rename(df, bene_id=BENE_ID, service_date=SRVC_DT, days_supply=DAYSSPLY, generic_name=GNN) %>%
    mutate(service_date=dmy(service_date)) %>%
    left_join(abx_class, by='generic_name') %>%
    select(-generic_name)
}

polish_year = function(y) {
  pde_in_fn = sprintf('../data/pde_%i.tsv', y)
  pde_out_fn = sprintf('pde_%i.tsv', y)

  read_tsv(pde_in_fn) %>%
    polish_pde() %>%
    write_tsv(pde_out_fn)

  bene_in_fn = sprintf('../data/bene_%i.tsv', y)
  bene_out_fn = sprintf('bene_%i.tsv', y)

  read_tsv(bene_in_fn, col_types=cols(zipcode='c')) %>%
    polish_bene() %>%
    write_tsv(bene_out_fn)
}

for (year in 2011:2014) polish_year(year)
