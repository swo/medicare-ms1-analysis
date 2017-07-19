# Proportions of diagnoses

cohort_ids = read_tsv('lm_cohort_ids.tsv') %$% bene_id
dx_class = read_tsv('../../data/fd_codes.tsv') %>%
  select(diagnosis=code, diagnosis_type)

dx = bind_rows(
  read_tsv('../../data/fd_dx_car_claims_2011.tsv') %>% mutate(time=0),
  read_tsv('../../data/fd_dx_op_claims_2011.tsv') %>% mutate(time=0),
  read_tsv('../../data/fd_dx_op_claims_2014.tsv') %>% mutate(time=1),
  read_tsv('../../data/fd_dx_car_claims_2014.tsv') %>% mutate(time=1)
) %>% rename(bene_id=BENE_ID) %>%
  filter(bene_id %in% cohort_ids) %>%
  mutate(from_date=dmy(from_date)) %>%
  left_join(dx_class, by='diagnosis')

baseline = count(dx, time, diagnosis_type) %>% ungroup() %>%
  arrange(time) %>% group_by(diagnosis_type) %>%
  summarize(dx_ratio=n[2]/n[1])

pde = bind_rows(
  read_tsv('../pde_2011.tsv') %>% mutate(time=0),
  read_tsv('../pde_2014.tsv') %>% mutate(time=1)
) %>% filter(bene_id %in% cohort_ids) %>%
  select(-days_supply) %>%
  mutate(pde_id=1:nrow(.))

# levofloxacin
window = 3

# age structure
p_age = bene %>% count(time, age) %>%
  group_by(time) %>% mutate(p_age=n/sum(n)) %>% ungroup() %>%
  select(time, age, p_age)

# unique days on which levo was dispensed to each bene
levo_days = pde %>% filter(antibiotic=='levofloxacin') %>%
  select(time, bene_id, service_date) %>%
  distinct() %>%
  mutate(day_id=1:nrow(.))

# same as above, but supplemented with a matching diagnosis
levo_days_dx = levo_days %>%
  left_join(dx, by=c('time', 'bene_id')) %>%
  filter(between(service_date - from_date, 0, window)) %>%
  group_by(day_id) %>% filter(from_date==max(from_date)) %>% ungroup() %>%
  select(day_id, diagnosis_type) %>%
  right_join(levo_days, by='day_id') %>%
  replace_na(list(diagnosis_type='none'))

# probability of a diagnosis (in a stretch of W days)
p_dxage = dx %>%
  select(time, bene_id, diagnosis_type, from_date) %>%
  group_by(time, bene_id, diagnosis_type) %>%
  summarize(n_dx_days=length(unique(from_date))) %>% ungroup() %>%
  mutate(p_dx=n_dx_days*window/365) %>%
  right_join(bene, by=c('time', 'bene_id')) %>%
  replace_na(list(p_dx=0)) %>%
  group_by(time, age, diagnosis_type) %>% summarize(p_dxage=mean(p_dx)) %>% ungroup()

p_rx = levo_days %>%
  count(time, bene_id) %>% ungroup() %>%
  right_join(bene, by=c('time', 'bene_id')) %>%
  replace_na(list(n=0)) %>%
  group_by(time, age) %>% summarize(p_rx=mean(n)/365) %>%
  ungroup()

matched_levo_days = dx %>% filter(diagnosis_type=='pneumonia') %>%
  left_join(levo_days, by=c('time', 'bene_id')) %>%
  filter(between(service_date - from_date, 0, window))

p_rxdxage = levo_days_dx %>%
  count(time, bene_id, diagnosis_type) %>% ungroup() %>%
  left_join(select(bene, time, bene_id, age), by=c('time', 'bene_id')) %>%
  group_by(time, age, diagnosis_type) %>%
  summarize(

p_rxdxage = levo_days %>% mutate(matched=day_id %in% matched_levo_days$day_id) %>%
  left_join(select(bene, time, bene_id, age), by=c('time', 'bene_id')) %>%
  count(time, age, matched) %>%
  group_by(time, age) %>% mutate(p_rxdxage=n/sum(n)) %>% ungroup() %>%
  filter(matched) %>%
  select(time, age, p_rxdxage)

p2 = p_rxdxage %>%
  left_join(p_dxage, by=c('time', 'age')) %>%
  left_join(p_age, by=c('time', 'age')) %>%
  select(time, age, p_rxdxage, p_dxage, p_age) %>%
  mutate(p=p_rxdxage*p_dxage*p_age)

levo_dx = levo_pde %>%
  left_join(dx, by=c('time', 'bene_id')) %>%
  filter(between(service_date - from_date, 0, window)) %>%
  # take the last (most recent) dx
  group_by(pde_id) %>% arrange(from_date) %>%
  mutate(dx_id=1:length(from_date)) %>% filter(dx_id==max(dx_id)) %>%
  select(-dx_id) %>% ungroup() %>%
  select(pde_id, diagnosis_type)
  
  
proportion_f = function(abx) {
  pde %>% filter(antibiotic==abx) %>%
    left_join(dx, by=c('bene_id', 'time')) %>%
    filter(between(service_date - from_date, 0, 3)) %>%
    count(time, diagnosis_type) %>% ungroup() %>%
    arrange(time) %>% group_by(diagnosis_type) %>%
    summarize(rx_ratio=n[2]/n[1]) %>%
    left_join(baseline, by='diagnosis_type') %>%
    mutate(p_ratio=rx_ratio/dx_ratio)
}

top_abx = c('azithromycin', 'levofloxacin',
            'amoxicillin/clavulanate', 'cephalexin', 'amoxicillin',
            'trimethoprim/sulfamethoxazole', 'clindamycin',
            'doxycycline', 'nitrofurantoin', 'ciprofloxacin')

res = lapply(top_abx, function(x) {
    proportion_f(x) %>%
      select(diagnosis_type, value=p_ratio) %>%
      mutate(abx=x)
  }) %>% bind_rows()

res %>%
  mutate(value=round(value, 2)) %>%
  spread(diagnosis_type, value) %>%
  kable