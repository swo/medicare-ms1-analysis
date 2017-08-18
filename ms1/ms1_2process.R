# Do all the heavy-lifting analyses so that the Rmd can just show the results
library(forcats)

# Load the census regions. Code them as factors so that Northeast is taken as
# the baseline in the linear models.
# NB: I put DC into the South
regions = read_tsv('../../db/census-regions/census-regions.tsv') %>%
  select(state, region) %>%
  mutate(region=factor(region, levels=c('Northeast', 'West', 'Midwest', 'South')))

dx_codes = read_tsv('../../data/fd_codes.tsv') %>%
  select(code, diagnosis_type)

# load data
subtract_min = function(x) x - min(x)

bene = read_tsv('data/bene.tsv') %>%
  mutate(heart_disease=AMI | ATRIALFB | CHF | ISCHMCHT) %>%
  mutate(is_female=sex=='female', is_white=race=='white', is_dual=buyin_months>0) %>%
  mutate(age=age-1) %>%
  left_join(regions, by='state') %>%
  mutate(year0=subtract_min(year),
         age0=subtract_min(age))

pde = read_tsv('data/pde.tsv')

dx_raw = read_tsv('data/dx.tsv') %>%
  left_join(dx_codes, by='code')

output_table = function(x, base) write_tsv(x, sprintf('tables/%s.tsv', base))
sem = function(x) sd(x) / sqrt(length(x))

summarize_totals = function(df) {
  summarize(df, n_bene=n(),
            mean_age=mean(age),
            sd_age=sd(age),
            sem_age=sem(age),
            mean_n_cc=mean(n_cc),
            sd_n_cc=sd(n_cc),
            sem_n_cc=sem(n_cc),
            n_female=sum(is_female),
            n_white=sum(is_white),
            n_dual=sum(is_dual),
            n_region_south=sum(region=='South'),
            n_region_west=sum(region=='West'),
            n_region_northeast=sum(region=='Northeast'),
            n_region_midwest=sum(region=='Midwest'))
}

# summary characteristics
totals = bene %>%
  group_by(year) %>%
  summarize_totals() %T>%
  output_table('totals')

n_unique_bene = bene$bene_id %>%
  unique %>% length %T>%
  write('tables/n_unique_bene.tsv')

# the top few inidividual antibiotics
claims_by_abx = pde %>%
  count(year, antibiotic) %>% ungroup() %>%
  left_join(select(totals, year, n_bene), by='year') %>%
  mutate(cpkp=n*1000/n_bene) %T>%
  output_table('claims_by_abx')

# keep track of which are the top 10 abx
top_abx = pde %>%
  count(antibiotic) %>%
  arrange(desc(n)) %$%
  head(antibiotic, 10)

# from here on, the denominators are taken from the beneficiary data grouped
# in the same way as the consumption data, so we can use the summarize_by function
summarize_by = function(bene, by) {
  group_by_(bene, 'year', by) %>%
    summarize(n_bene=n(), n_claims=sum(n_claims)) %>%
    ungroup() %>%
    mutate(cpkp=n_claims*1000/n_bene) %>%
    output_table(paste0('claims_by_', by))
}

summarize_by(bene, 'sex')
summarize_by(bene, 'region')
summarize_by(bene, 'race')
summarize_by(bene, 'age')

# run a model and tidy
glm_f = function(df, frmla, ...) eval(substitute(function(df2) glm(frmla, data=df2, ...)))(df)
model_f = function(df, frmla, name) {
  glm_f(df, frmla, family='poisson') %>%
    tidy %>%
    mutate(name=name)
}

abx_frmla = y ~ year0 + age0 + n_cc + is_female + is_white + is_dual + region

# get models and predictions for all drugs combined
overall_model = bene %>%
  rename(y=n_claims) %>%
  model_f(abx_frmla, 'overall')

# model overall consumption, but look at individual populations
covariates = c('age0', 'is_female', 'race', 'is_dual', 'region')
reduced_model = function(df, missing_term, name) {
  remaining_covariates = covariates[-which(covariates == missing_term)] %>% str_c(collapse=' + ')
  frmla = str_interp("y ~ year0 + n_cc + ${remaining_covariates}") %>% as.formula
  
  df %>%
    rename(y=n_claims) %>%
    model_f(frmla, name)
}

bind_rows(
  bene %>% filter(between(age, 65, 75)) %>% reduced_model('age0', 'age65_75'),
  bene %>% filter(between(age, 76, 85)) %>% reduced_model('age0', 'age76_85'),
  bene %>% filter(between(age, 86, 95)) %>% reduced_model('age0', 'age86_95'),
  bene %>% filter(is_female) %>% reduced_model('is_female', 'female'),
  bene %>% filter(!is_female) %>% reduced_model('is_female', 'male'),
  bene %>% filter(race=='white') %>% reduced_model('race', 'white'),
  bene %>% filter(race=='black') %>% reduced_model('race', 'black'),
  bene %>% filter(race=='Hispanic') %>% reduced_model('race', 'Hispanic'),
  bene %>% filter(race=='other') %>% reduced_model('race', 'other'),
  bene %>% filter(region=='South') %>% reduced_model('region', 'South'),
  bene %>% filter(region=='Midwest') %>% reduced_model('region', 'Midwest'),
  bene %>% filter(region=='West') %>% reduced_model('region', 'West'),
  bene %>% filter(region=='Northeast') %>% reduced_model('region', 'Northeast')
) %>%
  output_table('model_overall_demography')

# get models for individual abx
single_f = function(abx, frmla) {
  pde %>%
    filter(antibiotic==abx) %>%
    count(year, bene_id) %>% ungroup() %>% rename(y=n) %>%
    right_join(bene, by=c('year', 'bene_id')) %>%
    replace_na(list(y=0)) %>%
    model_f(frmla, abx)
}

single_models = lapply(top_abx, function(a) single_f(a, abx_frmla))

# report the model parameters
bind_rows(overall_model, single_models) %>%
  output_table('model_abx')

# trends in all fluoroquinolones
model_fq = pde %>%
  filter(str_detect(antibiotic, 'floxacin')) %>%
  count(year, bene_id) %>% ungroup() %>% rename(y=n) %>%
  right_join(bene, by=c('year', 'bene_id')) %>% replace_na(list(y=0)) %>%
  model_f(abx_frmla, 'fluoroquinolone') %T>%
  output_table('model_fq')

# special heart disease models
# how many beneficiaries have heart disease?
bene %>%
  count(year, heart_disease) %>%
  output_table('counts_heartdisease')

# did azithro consumption decrease faster in benes with heart disease?
heart_frmla = y ~ year0*heart_disease + age0 + n_cc + is_female + is_white + is_dual + region
model_azithro = single_f('azithromycin', heart_frmla) %T>%
  output_table('model_azithro')

model_levo = single_f('levofloxacin', heart_frmla) %T>%
  output_table('model_levo')

## DIAGNOSES ##
# count diagnoses
dx_counts = dx_raw %>%
  select(dx_id, year, diagnosis_type) %>%
  distinct() %>%
  count(year, diagnosis_type) %T>%
  output_table('dx_counts')

# count delays
dx_delays = dx_raw %>%
  count(delay) %T>%
  output_table('dx_delay')

# for each diagnosis, how often is each drug used?
dx_freqs = dx_raw %>%
  count(year, diagnosis_type, antibiotic) %>% ungroup() %T>%
  output_table('dxrx_fractions')

# create a table of dx (rows) with multiple abx as columns
dxrx = dx_raw %>%
  select(year, bene_id, from_date, diagnosis_type, antibiotic) %>%
  mutate(antibiotic=fct_other(antibiotic, keep=c('none', top_abx))) %>%
  distinct() %>%
  mutate(dummy=TRUE) %>% spread(antibiotic, dummy, fill=FALSE) %>%
  left_join(bene, by=c('year', 'bene_id'))

dx_types = unique(dxrx$diagnosis_type)
dx_replace = rep(FALSE, length(dx_types)) %>% as.list %>% setNames(dx_types)

dxrx_binary = dxrx %>%
  select(year, bene_id, diagnosis_type) %>%
  distinct() %>%
  mutate(dummy=TRUE) %>%
  spread(diagnosis_type, dummy) %>%
  right_join(bene, by=c('year', 'bene_id')) %>%
  replace_na(dx_replace)

# risk ratio
rr = function(or, p) or / (1 - p + (p * or))

# adjusted OR for diagnoses
# i.e., logistic regression bene_had_dx_in_year? ~ year + covariates
dx_adjust_f = function(dxrx, dx) {
  frmla = as.formula(str_interp("${dx} ~ year0 + age0 + n_cc + is_white + is_dual + is_female + region"))
  m = glm(frmla, data=dxrx, family='binomial')
  yearly_risks = dxrx %>%
    mutate(risk=predict(m, type='response')) %>%
    group_by(year0) %>%
    summarize(risk=mean(risk)) %>%
    mutate(term=str_interp("yearly_risk_year0=${year0}"), estimate=risk) %>%
    select(term, estimate)
  
  bind_rows(tidy(m), yearly_risks) %>%
    mutate(diagnosis_type=dx)
}

dx_res = lapply(dx_types, function(x) dx_adjust_f(dxrx_binary, x)) %>%
  bind_rows() %T>%
  output_table('model_dx')

# for a diagnosis type dx and an antibiotic, what is the trend in the probability
# that that drug will be used for that diagnosis?
# i.e., logistic regression drug_used? ~ dx_present?*year + covariates
# we also ask about "any diagnosis", i.e., drug_used? ~ year + covariates
dxrx_f = function(dxrx, dx, abx) {
  if (dx=='any_dx') {
    filter_f = function(df) df
  } else {
    filter_f = function(df) filter(df, diagnosis_type==dx)
  }
  
  frmla = as.formula(str_interp("`${abx}` ~ year0 + age0 + n_cc + is_white + is_dual + is_female + region"))
  
  dxrx %>%
    filter_f() %>%
    glm_f(frmla, family='binomial') %>%
    tidy %>%
    mutate(diagnosis_type=dx, antibiotic=abx)
}

dxrx_res = crossing(diagnosis_type=dx_types,
                    antibiotic=c('azithromycin', 'levofloxacin', 'amoxicillin/clavulanate', 'cephalexin', 'ciprofloxacin')) %>%
  rowwise() %>%
  do(dxrx_f(dxrx, .$diagnosis_type, .$antibiotic)) %>%
  ungroup() %T>%
  output_table('model_dxrx')
