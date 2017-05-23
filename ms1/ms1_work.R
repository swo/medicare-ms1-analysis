# Do all the heavy-lifting analyses so that the Rmd can just show the results

import::from(lazyeval, interp)
import::from(broom, tidy)
import::from(stringr, str_replace)

# Load the census regions. Code them as factors so that Northeast is taken as
# the baseline in the linear models.
# NB: I put DC into the South
regions = read_tsv('../../db/census-regions/census-regions.tsv') %>%
  select(state, region)

condition_names = read_tsv('../chronic-conditions/names.txt',
                           col_names=c('condition', 'condition_name'))

load_data = function(year) {
  bene = read_tsv(sprintf('../bene_%i.tsv', year)) %>%
    mutate(is_female=sex=='female', is_white=race=='white', is_dual=buyin_months>0) %>%
    left_join(regions, by='state') %>%
    mutate(year=year)
  
  pde = read_tsv(sprintf('../pde_%i.tsv', year)) %>%
    select(bene_id, antibiotic, antibiotic_class) %>%
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
sem = function(x) sd(x) / length(x)

# summary characteristics
totals = bene %>%
  group_by(year) %>%
  summarize(n_bene=n(),
            mean_age=mean(age),
            sem_age=sem(age),
            mean_n_cc=mean(n_cc),
            sem_n_cc=sem(n_cc),
            n_female=sum(is_female),
            n_white=sum(is_white),
            n_dual=sum(is_dual),
            n_region_south=sum(region=='South'),
            n_region_west=sum(region=='West'),
            n_region_northeast=sum(region=='Northeast'),
            n_region_midwest=sum(region=='Midwest')) %T>% 
  output_table('totals')
            
# the top few inidividual antibiotics
pde %>%
  count(year, antibiotic) %>%
  ungroup() %>%
  left_join(select(totals, year, n_bene), by='year') %>%
  mutate(cpkp=n*1000/n_bene) %T>%
  output_table('claims_by_abx')

# antibiotic classes
pde %>%
  count(year, antibiotic_class) %>%
  ungroup() %>%
  left_join(select(totals, year, n_bene), by='year') %>%
  mutate(cpkp=n*1000/n_bene) %T>%
  output_table('claims_by_class')

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
  summarize_by('race') %T>%
  output_table('claims_by_race')

# then by age
# (the input to summarize_by is modified, since we want to group by age *group*,
# not just by raw age)
claims_by_age = bene %>%
  mutate(age_group=case_when(
    between(.$age, 66, 76) ~ '66-76',
    between(.$age, 77, 86) ~ '77-86',
    between(.$age, 87, 96) ~ '87-96')) %>%
  summarize_by('age_group') %>%
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
#single_condition_models = lapply(conditions, function(cnd) single_condition_model(cnd)) %>%
#  bind_rows() %T>%
#  output_table('single_condition_models')

filter(bene, age >= 67) %>%
  lm(n_claims ~ year, data=.) %>%
  tidy %T>%
  output_table('model_year_only')

filter(bene, age >= 67) %>%
  lm(n_claims ~ year + age + is_female + is_dual + is_white + n_cc, data=.) %>%
  tidy %T>%
  output_table('model_plus')

#all_condition_formula = formula(paste0('n_claims ~ year + age + is_female + is_dual + is_white + region + ', paste0(conditions, collapse=' + ')))

# do a linear model to guess the coefficients
#all_condition_model = filter(bene, age >= 67) %>%
#  lm(all_condition_formula, data=.) %>%
#  tidy %T>%
#  output_table('model_all_condition')