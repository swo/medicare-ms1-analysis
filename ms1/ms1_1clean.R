# Process the bene, pde, and dx data
# so it's ready for analysis

# Keep chronic conditions that have a single reference year. Exclude Alzheimer's,
# which is a subset of Alz. & dementia
condition_names = read_tsv('../chronic-conditions/names.tsv') %>% filter(condition_code != 'ALZH')
names_1yr_ccs = condition_names %>% filter(condition_ref_years==1) %$% condition_code
sum_1yr_cc = paste0(names_1yr_ccs, collapse=' + ')

# window of days that you can look ahead from an rx to a dx
rxdx_window = 3

load_data = function(year) {
  # pde data
  pde = read_tsv(sprintf('../pde_%i.tsv', year)) %>%
    select(bene_id, service_date, antibiotic) %>%
    distinct() %>%
    mutate(year=year, pde_id=1:n())

  # total numbers of PDEs
  n_claims = count(pde, bene_id) %>% rename(n_claims=n)

  # count chronic conditions. 'n_cc' means one-year cc's
  cc = read_feather(sprintf('../cc_%i.feather', year)) %>%
    mutate_(n_cc=sum_1yr_cc) %>%
    ungroup()

  # join the cc, summary PDE, and dx data into bene
  bene = read_tsv(sprintf('../bene_%i.tsv', year)) %>%
    mutate(year=year) %>%
    left_join(cc, by='bene_id') %>%
    left_join(n_claims, by='bene_id') %>% replace_na(list(n_claims=0))
  
  # diagnosis claims
  dx = read_tsv(sprintf('../../data/tmp_hcpcs_%i.tsv', year)) %>%
    select(bene_id=BENE_ID, code=dx, from_date=dt1) %>%
    # filter for bene's we have
    semi_join(bene, by='bene_id') %>%
    mutate(from_date=dmy(from_date), year=year, dx_id=1:n())
  
  # associate a PDE with each diagnosis, if possible
  dx_pde = dx %>%
    left_join(pde, by=c('year', 'bene_id')) %>%
    mutate(delay=service_date - from_date) %>%
    filter(between(delay, 0, rxdx_window)) %>%
    select(dx_id, antibiotic, delay) %>%
    right_join(dx, by='dx_id') %>%
    replace_na(list(antibiotic='none'))

  # take the bene, pde (with diagnoses) and diagnoses (separately)
  list(bene=bene, pde=pde, dx=dx_pde)
}

dat = lapply(2011:2014, load_data)

save_dat = function(dat, name) {
  lapply(dat, function(x) getElement(x, name)) %>%
    bind_rows() %>%
    write_tsv(str_interp('data/${name}.tsv'))
}

lapply(c('bene', 'pde', 'dx'), function(x) save_dat(dat, x))
