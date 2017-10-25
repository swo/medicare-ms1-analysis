#!/usr/bin/env Rscript

# Process the bene, pde, and dx data so it's ready for analysis

# lower the names of a tibble
lower_names = function(df) {
  ns = names(df)
  rename(df, !!!setNames(as.list(ns), sapply(ns, tolower)))
}

regions = read_tsv('db/census-regions/census-regions.tsv') %>%
  select(state, region)

load_data = function(year) {
  # pde data
  pde = read_tsv(sprintf('raw_data/pde_%i.tsv', year)) %>%
    lower_names() %>%
    mutate(pde_date=dmy(pde_date)) %>%
    select(bene_id, pde_id, pde_date, antibiotic, fill_num) %>%
    # for each bene, date, and drug, keep only one record, giving preference to first fills.
    group_by(bene_id, pde_date, antibiotic) %>%
    summarize(pde_id=min(pde_id), fill_num=min(fill_num)) %>%
    ungroup() %>%
    mutate(year=year) %>%
    select(year, bene_id, pde_id, antibiotic, fill_num)

  # total numbers of PDEs
  n_claims = count(pde, bene_id) %>% rename(n_claims=n)
  n_claims_firstfill = pde %>%
    filter(fill_num==0) %>%
    count(bene_id) %>%
    rename(n_claims_firstfill=n)

  # join the summary PDE
  bene = read_tsv(sprintf('raw_data/bene_%i.tsv', year)) %>%
    lower_names() %>%
    left_join(regions, by='state') %>%
    left_join(n_claims, by='bene_id') %>%
    left_join(n_claims_firstfill, by='bene_id') %>%
    replace_na(list(n_claims=0, n_claims_firstfill=0)) %>%
    mutate(year=year, dual=buyin_mo>0) %>%
    select(year, bene_id, age, sex, race, dual, n_cc, region, n_claims, n_claims_firstfill)

  # take the bene, pde (with diagnoses) and diagnoses (separately)
  list(bene=bene, pde=pde)
}

save_dat = function(dat, name) {
  lapply(dat, function(x) getElement(x, name)) %>%
    bind_rows() %>%
    write_tsv(str_interp('data/${name}.tsv'))
}

dat = lapply(2011:2014, load_data)
lapply(c('bene', 'pde'), function(x) save_dat(dat, x))
