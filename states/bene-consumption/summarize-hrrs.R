#!/usr/bin/env Rscript

summarize_year = function(year, n_max=Inf) {
  hrr_fn = sprintf('../../../db/hrr/hrr_%i.tsv', year)
  usage_fn = sprintf('../../bene_usage_%i.tsv', year)
  output_fn = sprintf('hrr_usage_%i.tsv', year)

  hrr = read_tsv(hrr_fn)

  usage = read_tsv(usage_fn, n_max=n_max) %>%
    filter(age >= 65) %>%
    left_join(hrr, by='zipcode')

  denom = usage %>%
    count(hrr) %>%
    rename(n_bene=n)

  usage %>%
    filter(!is.na(antibiotic)) %>%
    group_by(hrr, antibiotic) %>%
    summarize(n_claims=sum(n_claims), n_days=sum(days_supply)) %>%
    left_join(denom, by='hrr') %>%
    mutate(cpkp=n_claims*1000/n_bene, did=n_days*1000/(365*n_bene)) %>%
    write_tsv(output_fn)
}

for (y in 2011:2014) summarize_year(y)
