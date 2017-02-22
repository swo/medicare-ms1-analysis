#!/usr/bin/env Rscript

# get the common chronic condition data
all_chronic = read_tsv("../../../data/cc_all.tsv", col_types='cdddd') %>%
  #swo> hack for column name
  rename(bene_id=bene)

summarize_year = function(year, n_max=Inf) {
  bene_fn = sprintf('../../bene_%i.tsv', year)
  pde_fn = sprintf('../../pde_%i.tsv', year)
  usage_fn = sprintf('state_usage_%i.tsv', year)

  year_column_name = sprintf('y%i', year)
  chronic = all_chronic %>%
    select_('bene_id', year_column_name) %>%
    rename_(.dots=setNames(year_column_name, 'chronic'))

  bene = read_tsv(bene_fn, n_max=n_max) %>%
    filter(!is.na(state), plan_coverage_months==12, age >= 65) %>%
    left_join(chronic, by='bene_id')

  bene_by_state = bene %>% count(state) %>% rename(n_ppl=n)

  pde = read_tsv(pde_fn, n_max=n_max) %>%
    inner_join(bene %>% select(bene_id), by='bene_id') # keep only PDEs for beneficiaries we have

  usage = pde %>%
    left_join(bene, by='bene_id') %>%
    group_by(state, antibiotic) %>%
    summarize(n_claims=n(), n_days=sum(days_supply)) %>%
    left_join(bene_by_state, by='state') %>%
    mutate(cpkp=n_claims*1000/n_ppl, did=n_days*1000/(365*n_ppl))

  write_tsv(usage, usage_fn)
}

for (y in 2011:2014) summarize_year(y)
