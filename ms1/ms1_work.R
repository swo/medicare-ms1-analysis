import::from(lazyeval, interp)
import::from(broom, tidy)
import::from(stringr, str_replace)

regions = read_tsv('../../db/census-regions/census-regions.tsv') %>%
  select(state, region)
condition_names = read_tsv('../chronic-conditions/names.txt',
                           col_names=c('condition', 'condition_name'))

load_data = function(year) {
  bene = read_tsv(sprintf('../bene_%i.tsv', year)) %>%
    filter(between(age, 66, 96), hmo_months==0, sex %in% c('male', 'female')) %>%
    select(-hmo_months) %>%
    mutate(is_female=sex=='female', is_white=race=='white') %>%
    left_join(regions, by='state') %>%
    mutate(year=year)
  
  pde = read_tsv(sprintf('../pde_%i.tsv', year)) %>%
    select(-service_date, -days_supply, -fill_number) %>%
    filter(bene_id %in% bene$bene_id) %>%
    mutate(year=year)
  
  # get total numbers of PDEs
  n_claims = pde %>% count(bene_id) %>% rename(n_claims=n)
  bene %<>%
    left_join(n_claims, by='bene_id') %>%
    replace_na(list(n_claims=0))
  
  cc = read_feather(sprintf('../cc_%i.feather', year)) %>%
    select(-ALZH) %>%
    mutate(n_cc=rowSums(select_(., '-bene_id')))
  
  bene %<>% left_join(cc, by='bene_id')
    
  list(bene=bene, pde=pde)
}

dat = lapply(2011:2014, load_data)
bene = lapply(dat, function(df) df$bene) %>% bind_rows
pde = lapply(dat, function(df) df$pde) %>% bind_rows
rm(dat)

output_table = function(x, base) write_tsv(x, sprintf('tables/tbl_%s.tsv', base))

# total number of beneficiaries and PDEs
totals = left_join(
  bene %>% count(year) %>% rename(n_bene=n),
  pde %>% count(year) %>% rename(n_pde=n),
  by='year') %>%
  mutate(cpkp=n_pde*1000/n_bene) %T>%
  output_table('totals')

# claims per 1,000 people by abx class
claims_by_class = pde %>%
  count(year, antibiotic_class) %>%
  ungroup() %>%
  left_join(select(totals, year, n_bene), by='year') %>%
  mutate(cpkp=n*1000/n_bene) %T>%
  output_table('claims_by_class')

# the top few inidividual antibiotics
claims_by_abx = pde %>%
  count(year, antibiotic) %>%
  ungroup() %>%
  left_join(select(totals, year, n_bene), by='year') %>%
  mutate(cpkp=n*1000/n_bene) %T>%
  output_table('claims_by_abx')

# from here on, the denominators are taken from the beneficiary data grouped
# in the same way as the consumption data, so we can use the summarize_by function
summarize_by = function(bene, by) {
  group_by_(bene, 'year', by) %>%
    summarize(n_bene=n(), n_claims=sum(n_claims)) %>%
    ungroup() %>%
    mutate(cpkp=n_claims*1000/n_bene)
}

claims_by_sex = summarize_by(bene, 'sex') %T>%
  output_table('claims_by_sex')
claims_by_region = summarize_by(bene, 'region') %T>%
  output_table('claims_by_region')

claims_by_race = bene %>%
  mutate(race=if_else(race %in% c('white', 'black', 'Hispanic', 'Asian'), race, 'other')) %>%
  summarize_by('race') %T>%
  output_table('claims_by_race')

# then by age
# (the input to summarize_by is modified, since we want to group by age *group*,
# not just by raw age)
claims_by_age = bene %>%
  mutate(group='all') %>%
  bind_rows(filter(., n_cc==0) %>% mutate(group='no_cc')) %>%
  group_by(year, group, age) %>%
  summarize(n_bene=n(),
    mean_n_claims=mean(n_claims),
    sd_n_claims=sd(n_claims),
    mean_n_cc=mean(n_cc),
    sd_n_cc=sd(n_cc)) %>%
  mutate(sem_n_claims=sd_n_claims/sqrt(n_bene),
         sem_n_cc=sd_n_cc/sqrt(n_bene)) %T>%
  output_table('claims_by_age')
  
cc_prevalence = bene %>%
  select(year, AMI:HYPOTH) %>%
  group_by(year) %>%
  summarize_all(mean) %>%
  gather('condition', 'prevalence', AMI:HYPOTH) %T>%
  output_table('cc_prevalence')

single_condition_model = function(condition_code) {
  filter(bene, age >= 67) %>%
    rename_(.dots=setNames(condition_code, 'has_condition')) %>%
    glm(n_claims ~ age + is_female + has_condition + year,
        family=quasipoisson(link='identity'), data=.) %>%
    tidy() %>%
    mutate(condition=condition_code)
}

conditions = c("AMI", "ALZHDMTA", "ATRIALFB", "CATARACT", "CHRNKIDN", "COPD", "CHF", "DIABETES", "GLAUCOMA", "HIPFRAC", "ISCHMCHT", "DEPRESSN", "OSTEOPRS", "RA_OA", "STRKETIA", "CNCRBRST", "CNCRCLRC", "CNCRPRST", "CNCRLUNG", "CNCRENDM", "ANEMIA", "ASTHMA", "HYPERL", "HYPERP", "HYPERT", "HYPOTH")
single_condition_models = lapply(conditions, function(cnd) single_condition_model(cnd)) %>%
  bind_rows() %T>%
  output_table('single_condition_models')

all_condition_formula = formula(paste0('n_claims ~ age + is_female + year + ', paste0(conditions, collapse=' + ')))
start = rep(0.1, length(all.vars(all_condition_formula)))
all_condition_model = filter(bene, age >= 67) %>%
  glm(all_condition_formula, family=quasipoisson(link='identity'), start=start, data=.) %>%
  tidy %T>%
  output_table('all_condition_model')
