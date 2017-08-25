#!/usr/bin/env Rscript

abx_class = read_tsv('abx_class.tsv')

polish_year = function(y) {
  # bene
  bene_in_fn = sprintf('../data/bene_%i.tsv', y)
  bene_out_fn = sprintf('bene_%i.tsv', y)

  read_tsv(bene_in_fn, col_types=cols(zipcode='c')) %>%
    select(bene_id=BENE_ID, state, sex, race, zipcode, age=AGE, buyin_months=BUYIN_MO) %>%
    write_tsv(bene_out_fn)

  # pde
  pde_in_fn = sprintf('../data/pde_%i.tsv', y)
  pde_out_fn = sprintf('pde_%i.tsv', y)

  read_tsv(pde_in_fn) %>%
    select(bene_id=BENE_ID, pde_id=PDE_ID, pde_date=SRVC_DT, days_supply=DAYSSPLY, gcdf=GCDF, generic_name=GNN) %>%
    mutate(pde_date=dmy(pde_date)) %>%
    left_join(abx_class, by='generic_name') %>%
    select(-generic_name) %>%
    write_tsv(pde_out_fn)

  # diagnoses: carrier and outpatient
  car_in_fn = sprintf('../data/car_%i.tsv', y)
  op_in_fn = sprintf('../data/op_%i.tsv', y)
  dx_out_fn = sprintf('dx_%i.tsv', y)

  bind_rows(
    read_tsv(car_in_fn) %>% select(bene_id=BENE_ID, claim_id=CLM_ID, dx_date=EXPNSDT1, dx_code=LINE_ICD_DGNS_CD, dx_category=diagnosis_category),
    read_tsv(op_in_fn) %>% select(bene_id=BENE_ID, claim_id=CLM_ID, dx_date=THRU_DT, dx_code=PRNCPAL_DGNS_CD, dx_category=diagnosis_category)
  ) %>%
    mutate(dx_date=dmy(dx_date)) %>%
    distinct() %>%
    write_tsv(dx_out_fn)

  # chronic conditions
  cc_fn = sprintf('../data/cc_%i.tsv', y)
  cc_out_fn = sprintf('cc_%i.feather', y)
  read_tsv(cc_fn) %>%
    rename(bene_id=BENE_ID) %>%
    mutate_at(vars(AMI:HYPOTH), function(x) x==3) %>%
    write_feather(cc_out_fn)
}

library(parallel)
cluster = makeCluster(3, type='FORK')

parLapply(cluster, 2011:2014, polish_year)

stopCluster(cluster)
