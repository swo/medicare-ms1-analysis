#!/usr/bin/env Rscript

summarize_year = function(year) {
  hrr_fn = sprintf('../../../db/hrr/hrr_%i.tsv', year)
  bene_fn = paste0('../../bene_', year, '.tsv')
  pde_fn = paste0('../../pde_', year, '.tsv')

  hrr = read_tsv(hrr_fn)

  bene = read_tsv(bene_fn) %>%
    filter(age >= 65, hmo_months==0) %>%
    select(bene_id, zipcode) %>%
    left_join(hrr, by='zipcode')

  bene_by_hrr = bene %>% count(hrr) %>% rename(n_ppl=n)

  pde = read_tsv(pde_fn)

  usage = pde %>%
    left_join(bene, by='bene_id') %>%
    group_by(hrr, antibiotic) %>%
    summarize(n_claims=n(), n_days=sum(days_supply)) %>%
    left_join(bene_by_hrr, by='hrr') %>%
    mutate(cpkp=n_claims*1000/n_ppl, did=n_days*1000/(365*n_ppl)) %>%
    mutate(year=year)

  usage
}

lapply(2011:2014, summarize_year) %>%
  bind_rows() %>%
  write_tsv('hrr_usage.tsv')
