#!/usr/bin/env Rscript

library(assertive)

# create summary beneficiary x abx file

# get the common chronic condition data
all_chronic = read_tsv("../data/cc_all.tsv", col_types='cdddd') %>%
  #swo> hack for column name
  rename(bene_id=bene)

regions = read_tsv('../db/census-regions/census-regions.tsv') %>%
  select(state, region)

antibiotic_class = read_tsv('abx_class.tsv') %>%
  select(antibiotic, antibiotic_class) %>%
  distinct

assert_has_no_duplicates(antibiotic_class$antibiotic)

summarize_year = function(year, n_max=Inf) {
  bene_fn = sprintf('bene_%i.tsv', year)
  pde_fn = sprintf('pde_%i.tsv', year)
  usage_fn = sprintf('bene_usage_%i.tsv', year)

  year_column_name = sprintf('y%i', year)
  chronic = all_chronic %>%
    select_('bene_id', year_column_name) %>%
    rename_(.dots=setNames(year_column_name, 'chronic'))

  bene = read_tsv(bene_fn, n_max=n_max) %>%
    mutate_at(vars(age, plan_coverage_months), as.integer) %>%
    filter(!is.na(state), plan_coverage_months==12, age >= 65) %>%
    left_join(chronic, by='bene_id') %>%
    #swo> filter for only those beneficiaries that appear in the chronic data?
    left_join(regions, by='state')
  
  assert_has_no_duplicates(bene$bene_id)

  pde = read_tsv(pde_fn, n_max=n_max) %>%
    mutate_at(vars(days_supply), as.integer) %>%
    left_join(antibiotic_class, by='antibiotic')
  
  usage = bene %>%
    left_join(pde, by='bene_id') %>%
    mutate(n_claims=if_else(is.na(antibiotic), 0L, 1L),
           n_days=if_else(is.na(days_supply), 0L, days_supply))

  write_tsv(usage, usage_fn)
}

for (y in 2011:2014) summarize_year(y)
