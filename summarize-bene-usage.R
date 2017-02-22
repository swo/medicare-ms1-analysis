#!/usr/bin/env Rscript

# create summary beneficiary x abx file

# get the common chronic condition data
all_chronic = read_tsv("../data/cc_all.tsv", col_types='cdddd') %>%
  #swo> hack for column name
  rename(bene_id=bene)

regions = read_tsv('../db/census-regions/census-regions.tsv') %>%
  select(state, region)

summarize_year = function(year, n_max=Inf) {
  bene_fn = sprintf('bene_%i.tsv', year)
  pde_fn = sprintf('pde_%i.tsv', year)
  usage_fn = sprintf('bene_usage_%i.tsv', year)

  year_column_name = sprintf('y%i', year)
  chronic = all_chronic %>%
    select_('bene_id', year_column_name) %>%
    rename_(.dots=setNames(year_column_name, 'chronic'))

  bene = read_tsv(bene_fn, n_max=n_max) %>%
    filter(!is.na(state), plan_coverage_months==12, age >= 65) %>%
    left_join(chronic, by='bene_id') %>%
    #swo> filter for only those beneficiaries that appear in the chronic data?
    left_join(regions, by='state')

  pde = read_tsv(pde_fn, n_max=n_max) %>%
    inner_join(bene %>% select(bene_id), by='bene_id') # keep only PDEs for beneficiaries we have

  usage = pde %>%
    group_by(bene_id, antibiotic) %>%
    summarize(n_claims=n(), n_days=sum(days_supply)) %>%
    gather('consumption_metric', 'consumption', n_claims:n_days) %T>%
    { . %>% select(antibiotic) %>% distinct %>% write_tsv('tmp') } %>%
    spread(antibiotic, consumption, fill=0) %>%
    left_join(bene, by='bene_id')

  write_tsv(usage, usage_fn)
}

for (y in 2011:2014) summarize_year(y)
