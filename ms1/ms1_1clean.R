# Process the bene, pde, and dx data
# so it's ready for analysis

# Keep chronic conditions that have a single reference year. Exclude Alzheimer's,
# which is a subset of Alz. & dementia
condition_names = read_tsv('../chronic-conditions/names.tsv') %>% filter(condition_code != 'ALZH')
names_1yr_ccs = condition_names %>% filter(condition_ref_years==1) %$% condition_code
sum_1yr_cc = paste0(names_1yr_ccs, collapse=' + ')

regions = read_tsv('../../db/census-regions/census-regions.tsv') %>%
  select(state, region)

load_data = function(year) {
  # pde data
  pde = read_tsv(sprintf('../pde_%i.tsv', year)) %>%
    select(bene_id, pde_date, antibiotic, fill_num) %>%
    # for each bene, date, and drug, keep only one record. give preference to first fills.
    group_by(bene_id, pde_date, antibiotic) %>%
    summarize(fill_num=min(fill_num)) %>%
    ungroup() %>%
    mutate(year=year, pde_id=1:n())

  # total numbers of PDEs
  n_claims = count(pde, bene_id) %>% rename(n_claims=n)
  n_claims_firstfill = pde %>%
    filter(fill_num==0) %>%
    count(bene_id) %>%
    rename(n_claims_firstfill=n)

  # count chronic conditions. 'n_cc' means one-year cc's
  cc = read_feather(sprintf('../cc_%i.feather', year)) %>%
    mutate_(n_cc=sum_1yr_cc) %>%
    ungroup() %>%
    select(bene_id, n_cc)

  # join the cc, summary PDE, and dx data into bene
  bene = read_tsv(sprintf('../bene_%i.tsv', year)) %>%
    mutate(year=year) %>%
    left_join(regions, by='state') %>%
    left_join(cc, by='bene_id') %>%
    left_join(n_claims, by='bene_id') %>%
    left_join(n_claims_firstfill, by='bene_id') %>%
    replace_na(list(n_claims=0, n_claims_firstfill=0))

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
