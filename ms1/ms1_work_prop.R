# Proportions of diagnoses

# Restrict analysis to beneficiaries in the cohort
cohort_ids = read_tsv('lm_cohort_ids.tsv') %$% bene_id
dx_class = read_tsv('../../data/fd_codes.tsv') %>%
  select(diagnosis=code, diagnosis_type)

# Load diagnoses from years 2011 ("0") and 2014 ("1")
dx = bind_rows(
  read_tsv('../../data/fd_dx_car_claims_2011.tsv') %>% mutate(time=0),
  read_tsv('../../data/fd_dx_op_claims_2011.tsv') %>% mutate(time=0),
  read_tsv('../../data/fd_dx_op_claims_2014.tsv') %>% mutate(time=1),
  read_tsv('../../data/fd_dx_car_claims_2014.tsv') %>% mutate(time=1)
) %>% rename(bene_id=BENE_ID) %>%
  filter(bene_id %in% cohort_ids) %>%
  mutate(from_date=dmy(from_date)) %>%
  left_join(dx_class, by='diagnosis') %>%
  left_join(select(bene, time, bene_id, age), by=c('time', 'bene_id')) %>%
  select(time, bene_id, age, from_date, diagnosis, diagnosis_type)

# load PDEs from those years
pde = bind_rows(
  read_tsv('../pde_2011.tsv') %>% mutate(time=0),
  read_tsv('../pde_2014.tsv') %>% mutate(time=1)
) %>% filter(bene_id %in% cohort_ids) %>%
  select(-days_supply) %>%
  mutate(pde_id=1:nrow(.)) %>%
  left_join(select(bene, time, bene_id, age), by=c('time', 'bene_id'))

# a window 
window = 7

# age structure
p_age = bene %>% count(time, age) %>%
  group_by(time) %>% mutate(p_age=n/sum(n)) %>% ungroup() %>%
  mutate(key=sprintf('p_age_t%i', time)) %>%
  select(age, key, p_age) %>%
  spread(key, p_age)

# a "narrow" age structure
# it excludes the old and young ages that don't appear in both
# the 2011 and 2014 data, so we can adjust for age structure
p_age_narrow = p_age %>%
  select(age, p_age=p_age_t0) %>%
  filter(between(age, 69, 92)) %>%
  mutate(p_age=p_age/sum(p_age))

n_dxage = dx %>%
  count(time, age, diagnosis_type) %>% ungroup()

# same as above, but supplemented with a matching diagnosis
pde_dx = pde %>%
  select(pde_id, time, bene_id, age, service_date, antibiotic) %>%
  left_join(dx, by=c('time', 'bene_id')) %>%
  filter(between(service_date - from_date, 0, window)) %>%
  group_by(pde_id) %>% filter(from_date==max(from_date)) %>% ungroup() %>%
  select(pde_id, diagnosis_type) %>%
  right_join(pde, by='pde_id') %>%
  replace_na(list(diagnosis_type='none'))
  
n_rxdxage = pde_dx %>%
  count(time, age, antibiotic, diagnosis_type) %>% ungroup()

proportion_f = function(abx) {
  dx_counts = n_dxage %>%
    mutate(key=sprintf('n_dx_t%i', time)) %>%
    select(age, diagnosis_type, key, n) %>%
    spread(key, n)
  rx_counts = n_rxdxage %>%
    filter(antibiotic==abx) %>%
    mutate(key=sprintf('n_rx_t%i', time)) %>%
    select(age, diagnosis_type, key, n) %>%
    spread(key, n)
  
  overall_rx_counts = n_rxdxage %>%
    filter(antibiotic==abx) %>%
    group_by(time, age) %>% summarize(n=sum(n)) %>% ungroup() %>%
    mutate(key=sprintf('n_rx_t%i', time)) %>%
    select(age, key, n) %>%
    spread(key, n) %>% mutate(r=n_rx_t1/n_rx_t0) %>%
    left_join(p_age, by='age') %>%
    mutate(r=(n_rx_t1/n_rx_t0)/(p_age_t1/p_age_t0)) %>%
    select(age, r) %>%
    mutate(diagnosis_type='overall')
  
  inner_join(dx_counts, rx_counts, by=c('age', 'diagnosis_type')) %>%
    mutate(r=(n_rx_t1/n_rx_t0)/(n_dx_t1/n_dx_t0)) %>%
    select(age, diagnosis_type, r) %>%
    bind_rows(overall_rx_counts) %>%
    filter(between(age, 69, 92)) %>%
    left_join(p_age_narrow, by='age') %>%
    group_by(diagnosis_type) %>%
    summarize(r=sum(r*p_age))
}

top_abx = c('azithromycin', 'levofloxacin',
            'amoxicillin/clavulanate', 'cephalexin', 'amoxicillin',
            'trimethoprim/sulfamethoxazole', 'clindamycin',
            'doxycycline', 'nitrofurantoin', 'ciprofloxacin')

res = lapply(top_abx, function(x) {
    proportion_f(x) %>%
      select(diagnosis_type, value=r) %>%
      mutate(abx=x)
  }) %>% bind_rows()

f_rxdx = n_rxdxage %>%
  group_by(antibiotic, diagnosis_type) %>%
  summarize(n=sum(n)) %>%
  group_by(antibiotic) %>%
  mutate(f_rxdx=n/sum(n)) %>%
  ungroup() %>%
  select(abx=antibiotic, diagnosis_type, f_rxdx)

res2 = res %>%
  left_join(f_rxdx, by=c('abx', 'diagnosis_type')) %>%
  mutate(f_rxdx=if_else(diagnosis_type=='overall', 1.0, f_rxdx)) %>%
  filter(!is.na(value)) %>%
  select(abx, diagnosis_type, p_rxdx_ratio=value, f_rxdx)

write_tsv(res2, 'tables/tbl_rxdx_ratios.tsv')

res_f = function(x) {
  res2 %>%
    filter(abx==x) %>% select(-abx) %>%
    arrange(desc(f_rxdx)) %>%
    mutate_if(is.numeric, function(x) round(x, 2)) %>%
    kable
}
