# Process the bene, pde, and dx data
# so it's ready for analysis

# Keep chronic conditions that have a single reference year. Exclude Alzheimer's,
# which is a subset of Alz. & dementia
condition_names = read_tsv('../chronic-conditions/names.tsv') %>% filter(condition_code != 'ALZH')
names_1yr_ccs = condition_names %>% filter(condition_ref_years==1) %$% condition_code
sum_1yr_cc = paste0(names_1yr_ccs, collapse=' + ')

# window of days that you can look ahead from an rx to a dx
rxdx_window = 3

# diagnosis types
fd = read_tsv('../../db/fd_icd/fd_categories.tsv') %>%
  mutate(dc_id=1:n())

load_data = function(year) {
  # pde data
  pde = read_tsv(sprintf('../pde_%i.tsv', year)) %>%
    select(bene_id, service_date, antibiotic) %>%
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
    left_join(n_claims, by='bene_id') %>%
    replace_na(list(n_claims=0))
  
  # diagnosis claims
  dx = read_tsv(str_interp("../../data/dx_${year}.tsv")) %>%
    select(bene_id=BENE_ID, code=dx_code, diagnosis_category, from_date=dt1) %>%
    # filter for bene's we have
    semi_join(bene, by='bene_id') %>%
    mutate(from_date=dmy(from_date), year=year) %>%
    distinct()
  
  # associate a PDE with each diagnosis, if possible
  pde_dx = pde %>%
    # keep all PDEs for which there are matching dx's; keep combos
    inner_join(dx, by=c('year', 'bene_id')) %>%
    # keep only where the delay is in a certain window
    mutate(delay=service_date - from_date) %>%
    filter(between(delay, 0, rxdx_window)) %>%
    # if there are multiple window-valid dx's per PDE,
    # pick the ones closest in time. if still tied,
    # pick the one with the highest-rated diagnosis category
    left_join(fd, by='diagnosis_category') %>%
    group_by(pde_id) %>%
    filter(delay==min(delay)) %>%
    filter(dc_id==min(dc_id)) %>%
    # and then just take the first one
    filter(row_number()==1) %>%
    ungroup() %>%
    select(pde_id, diagnosis_category, code, delay) %>%
    right_join(pde) %>%
    replace_na(list(diagnosis_category='none'))

  # take the bene, pde (with diagnoses) and diagnoses (separately)
  list(bene=bene, pde=pde_dx, dx=dx)
}

dat = lapply(2011:2014, load_data)

save_dat = function(dat, name) {
  lapply(dat, function(x) getElement(x, name)) %>%
    bind_rows() %>%
    write_tsv(str_interp('data/${name}.tsv'))
}

lapply(c('bene', 'pde', 'dx'), function(x) save_dat(dat, x))
