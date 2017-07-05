sum_over = function(target_columns) {
  dummy_vars = sapply(1:length(target_columns), function(i) sprintf('x%i', i))
  dummy_vars_string = paste0(dummy_vars, collapse=',')
  dummy_vars_list = lapply(target_columns, as.name) %>% setNames(dummy_vars)
  f = as.formula(paste0(c('~sum(', dummy_vars_string, ')'), collapse=''))
  lazyeval::interp(f, .values=dummy_vars_list)
}

regions = read_tsv('../../db/census-regions/census-regions.tsv') %>%
  select(state, region)

condition_names = read_tsv('../chronic-conditions/names.tsv') %>% filter(condition_code != 'ALZH')
names_1yr_ccs = condition_names %>% filter(condition_ref_years==1) %$% condition_code

dx_codes = read_tsv('../sex-difference/sex_codes.tsv')

load_data = function(year) {
  # pde data
  pde = read_tsv(sprintf('../pde_%i.tsv', year)) %>%
    select(bene_id, antibiotic, days_supply) %>%
    group_by(bene_id, antibiotic) %>%
    summarize(n_claims=n(), days_supply=sum(days_supply)) %>%
    ungroup() %>%
    mutate(year=year)

  # diagnosis claims
  dx_car_fn = sprintf('../../data/dx_car_claims_%i.tsv', year)
  dx_op_fn = sprintf('../../data/dx_op_claims_%i.tsv', year)

  dx = bind_rows(read_tsv(dx_car_fn),
                 read_tsv(dx_op_fn)) %>%
    mutate(date=dmy(from_date)) %>%
    rename(bene_id=BENE_ID, code=diagnosis) %>%
    left_join(dx_codes, by='code')

  dx_presence = dx %>%
    group_by(bene_id) %>%
    summarize(rc='acute_rc' %in% diagnosis_type,
              uti='uti' %in% diagnosis_type)

  count_car_fn = sprintf('../../data/count_car_claims_%i.tsv', year)
  count_op_fn = sprintf('../../data/count_op_claims_%i.tsv', year)
  n_op_days = bind_rows(read_tsv(count_car_fn), read_tsv(count_op_fn)) %>%
    rename(bene_id=BENE_ID) %>%
    group_by(bene_id) %>%
    summarize(n_op_days=sum(n_from_date))

  # total numbers of PDEs
  n_claims = pde %>% group_by(bene_id) %>% summarize(n_claims=sum(n_claims))

  # count chronic conditions. 'n_cc' means one-year cc's
  cc = read_feather(sprintf('../cc_%i.feather', year)) %>%
    rowwise() %>%
    mutate_(n_cc=sum_over(names_1yr_ccs)) %>%
    ungroup()

  # join the cc, summary PDE, and dx data into bene
  bene = read_tsv(sprintf('../bene_%i.tsv', year)) %>%
    mutate(is_female=sex=='female', is_white=race=='white', is_dual=buyin_months>0) %>%
    left_join(regions, by='state') %>%
    mutate(year=year) %>%
    left_join(cc, by='bene_id') %>%
    left_join(n_op_days, by='bene_id') %>% replace_na(list(n_op_days=0)) %>%
    left_join(n_claims, by='bene_id') %>% replace_na(list(n_claims=0)) %>%
    left_join(dx_presence, by='bene_id') %>% replace_na(list(rc=FALSE, uti=FALSE))

  list(bene=bene, pde=pde, dx=dx)
}

years = 2011:2014
n_years = length(years)
dat = lapply(years, load_data)

bene = lapply(dat, function(df) df$bene) %>% bind_rows %>%
  group_by(bene_id) %>%
  mutate(n_years=n(),
         start_age=min(age),
         in_cohort=n_years==4) %>%
  ungroup() %T>%
  write_tsv('data/bene.tsv')

pde = lapply(dat, function(df) df$pde) %>% bind_rows %>%
  left_join(distinct(select(bene, bene_id, in_cohort)), by='bene_id') %T>%
  write_tsv('data/pde.tsv')

dx = lapply(dat, function(df) df$dx) %>% bind_rows() %>%
  left_join(distinct(select(bene, bene_id, in_cohort)), by='bene_id') %T>%
  write_tsv('data/dx.tsv')

# save diagnostic counts
dx %>%
  filter(in_cohort) %>%
  count(code) %>%
  write_tsv('data/tbl_dx_counts.tsv')

# save top antibiotics
all_abx = pde %>%
  count(antibiotic) %T>%
  write_tsv('data/tbl_top_abx.tsv')

top_abx = all_abx %>% arrange(desc(n)) %$% head(antibiotic, 10)

# consumption of drugs by sex
usage_by_sex = function(drug) {
  bene %>%
    filter(in_cohort) %>%
    left_join(pde %>% filter(antibiotic==drug) %>% count(bene_id), by='bene_id') %>%
    replace_na(list(n=0)) %>%
    group_by(sex) %>%
    summarize(n_bene=n(),
              n_using=sum(n>0),
              n_claims=sum(n),
              mean_n_claims=mean(n),
              sd_n_claims=sd(n)) %>%
    mutate(se_n_claims=sd_n_claims/sqrt(n_bene),
           drug=drug)
}

lapply(top_abx, usage_by_sex) %>% bind_rows %>%
  write_tsv('data/tbl_usage_by_sex.tsv')

# healthcare_usage_by_sex
dx %>%
  count(bene_id, diagnosis_type) %>%
  ungroup() %>%
  right_join(bene, by='bene_id') %>%
  replace_na(list(n=0)) %>%
  group_by(sex) %>%
  summarize(total_n_days=sum(n),
            mean_n_days=mean(n),
            n_using=sum(n>0),
            n_bene=n()) %>%
  write_tsv('data/tbl_hc_by_sex.tsv')

cipro_dat = pde %>%
  filter(in_cohort, antibiotic=='ciprofloxacin') %>%
  select(bene_id, year, cipro=n_claims) %>%
  right_join(filter(bene, in_cohort), by=c('bene_id', 'year')) %>%
  replace_na(list(cipro=0)) %>%
  group_by(bene_id) %>%
  summarize(is_female=any(is_female),
            uti=any(uti),
            cipro=sum(cipro))

model_f = function(df, f) glm(f, data=df, family='poisson') %>% tidy

cipro_results = bind_rows(model_f(cipro_dat, cipro ~ is_female) %>% mutate(model='base'),
                          model_f(cipro_dat, cipro ~ is_female + uti) %>% mutate(model='expl'))
write_tsv(cipro_results, 'data/tbl_model_cipro.tsv')

bactrim_dat = pde %>%
  filter(in_cohort, antibiotic=='trimethoprim/sulfamethoxazole') %>%
  select(bene_id, year, bactrim=n_claims) %>%
  right_join(filter(bene, in_cohort), by=c('bene_id', 'year')) %>%
  replace_na(list(bactrim=0)) %>%
  group_by(bene_id) %>%
  summarize(is_female=any(is_female),
            uti=any(uti),
            bactrim=sum(bactrim))

bactrim_results = bind_rows(model_f(bactrim_dat, bactrim ~ is_female) %>% mutate(model='base'),
                            model_f(bactrim_dat, bactrim ~ is_female + uti) %>% mutate(model='expl'))
write_tsv(bactrim_results, 'data/tbl_model_tmpsmx.tsv')

nitro_dat = pde %>%
  filter(in_cohort, antibiotic=='nitrofurantoin') %>%
  select(bene_id, year, n=n_claims) %>%
  right_join(filter(bene, in_cohort), by=c('bene_id', 'year')) %>%
  replace_na(list(n=0)) %>%
  group_by(bene_id) %>%
  summarize(is_female=any(is_female),
            uti=any(uti),
            n=sum(n))

nitro_results = bind_rows(model_f(nitro_dat, n ~ is_female) %>% mutate(model='base'),
                          model_f(nitro_dat, n ~ is_female + uti) %>% mutate(model='expl'))
write_tsv(nitro_results, 'data/tbl_model_nitro.tsv')

azithro_dat = pde %>%
  filter(in_cohort, antibiotic=='azithromycin') %>%
  select(bene_id, year, n=n_claims) %>%
  right_join(filter(bene, in_cohort), by=c('bene_id', 'year')) %>%
  replace_na(list(n=0)) %>%
  group_by(bene_id) %>%
  summarize(is_female=any(is_female),
            rc=any(rc),
            copd=any(COPD),
            n=sum(n))

azithro_results = bind_rows(model_f(azithro_dat, n ~ is_female) %>% mutate(model='base'),
                            model_f(azithro_dat, n ~ is_female + rc + copd) %>% mutate(model='expl'))

write_tsv(azithro_results, 'data/tbl_model_azithro.tsv')

# Consumption by age and sex
bene %>%
  filter(in_cohort) %>%
  group_by(bene_id) %>%
  summarize(is_female=any(is_female),
            n=sum(n_claims),
            age=min(age)) %>%
  write_tsv('data/tbl_consumption_age_sex.tsv')

# Sex by age
bene %>% filter(in_cohort) %>%
  count(age, sex) %>%
  write_tsv('data/tbl_sex_by_age.tsv')

# Margins of azithromycin consumption
pde %>%
  filter(antibiotic=='azithromycin') %>%
  group_by(bene_id) %>% summarize(azi=sum(n_claims)) %>%
  right_join(bene %>% filter(in_cohort) %>% select(bene_id, sex) %>% distinct,
    by='bene_id') %>% replace_na(list(azi=0)) %>%
  count(sex, azi) %>%
  group_by(sex) %>%
  mutate(f_ppl=n/sum(n), cum_f_ppl=cumsum(f_ppl),
         cum_cons=cumsum(azi*n), cum_f_cons=cum_cons/max(cum_cons)) %>%
  write_tsv('data/tbl_azithro_evenness.tsv')
