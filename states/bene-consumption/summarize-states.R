#!/usr/bin/env Rscript

# get the common chronic condition data
#all_chronic = read_tsv("../../../data/cc_all.tsv", col_types='cdddd') %>%
  #swo> hack for column name
  #rename(bene_id=bene)

summarize_year = function(year) {
  bene_fn = sprintf('../../bene_%i.tsv', year)
  pde_fn = sprintf('../../pde_%i.tsv', year)
  #usage_fn = sprintf('state_usage_%i.tsv', year)

  #year_column_name = sprintf('y%i', year)
  #chronic = all_chronic %>%
    #select_('bene_id', year_column_name) %>%
    #rename_(.dots=setNames(year_column_name, 'chronic'))

  bene = read_tsv(bene_fn) %>%
    filter(hmo_months==0, age >= 65) %>%
    select(bene_id, state)

  bene_by_state = bene %>% count(state) %>% rename(n_ppl=n)

  pde = read_tsv(pde_fn)

  usage = pde %>%
    left_join(bene, by='bene_id') %>%
    group_by(state, antibiotic) %>%
    summarize(n_claims=n(), n_days=sum(days_supply)) %>%
    left_join(bene_by_state, by='state') %>%
    mutate(cpkp=n_claims*1000/n_ppl, did=n_days*1000/(365*n_ppl)) %>%
    mutate(year=year)

  usage
}

lapply(2011:2014, summarize_year) %>%
  bind_rows() %>%
  write_tsv('state_usage.tsv')
