#!/usr/bin/env Rscript

# load data
# NB: I put DC into the South

output_table = function(x, base) write_tsv(x, sprintf('tables/%s.tsv', base))

# bene
bene = read_tsv('data/bene.tsv') %>%
  mutate(region=factor(region, levels=c('Northeast', 'West', 'Midwest', 'South')),
         sex=factor(sex, levels=c('male', 'female')),
         race=factor(race, levels=c('white', 'black', 'Hispanic', 'other'))) %>%
  mutate(age=age-1)

read_years = function(years, template) {
  lapply(years, function(y) {
    read_tsv(sprintf(template, y)) %>%
      mutate(year=y)
  }) %>%
    bind_rows()
}

glm_f = function(df, frmla, ...) eval(substitute(function(df2) glm(frmla, data=df2, ...)))(df)
sandwich_tidy = function(m) tidy(lmtest::coeftest(m, vcov=sandwich::vcovHC(m, type='HC3')))
  
frmla = y ~ year + age + n_cc + sex + race + dual + region

# count total dx's and hcpcs counts
read_years(2011:2014, '../../data/dx_hcpcs_count_%i.tsv') %>%
  output_table('hcpcs_count')

# count dx by categories
read_years(2011:2014, '../../data/dx_cat_count_%i.tsv') %>%
  rename(dx_cat=diagnosis_category) %>%
  left_join(bene %>% count(year) %>% rename(n_bene=n), by='year') %>%
  output_table('dx_cat_counts')

# count many PDEs have any dx or E&M dx upstream
read_years(2011:2014, '../../data/dx_any_pde_%i.tsv') %>%
  mutate_at(vars(has_dx, has_em_dx), as.logical) %>%
  count(year, has_dx, has_em_dx) %>%
  output_table('pde_with_dx_upstream')

# which dx's contribute to each abx?
dx_from_pde = read_years(2011:2014, '../../data/dx_from_pde_%i.tsv') %>%
  rename(bene_id=BENE_ID, pde_id=PDE_ID, dx_cat=diagnosis_category)

dxs = unique(dx_from_pde$dx_cat)

top_abx = c('azithromycin', 'ciprofloxacin', 'amoxicillin', 'cephalexin',
            'trimethoprim/sulfamethoxazole', 'levofloxacin', 'amoxicillin/clavulanate',
            'doxycycline', 'nitrofurantoin', 'clindamycin')

# get a denominator for total number of PDEs
pde_denom_2011 = dx_from_pde %>%
  filter(year==2011) %>%
  filter(antibiotic %in% top_abx) %>%
  select(pde_id, antibiotic) %>%
  distinct()

# compute the table just from 2011
dx_from_pde_table = crossing(antibiotic=top_abx, dx_cat=dxs) %>%
  group_by(antibiotic, dx_cat) %>%
  do((function(a, dx) {
    # how what fraction of PDEs for this abx had this dx upstream?
    pde_denom_2011 %>%
      filter(antibiotic==a) %>%
      left_join(filter(dx_from_pde, antibiotic==a, dx_cat==dx), by='pde_id') %>%
      mutate(present=!is.na(dx_cat)) %>%
      count(present)
  })(.$antibiotic, .$dx_cat)) %T>%
  output_table('dx_from_pde')

# appropriateness
fd = read_tsv('../../db/fd_icd/fd_categories.tsv') %>%
  select(dx_cat=diagnosis_category, tier)

pde_approp = dx_from_pde %>%
  left_join(fd, by='dx_cat') %>%
  # "Remaining codes not listed elsewhere" are tier 3
  replace_na(list(tier=3)) %>%
  group_by(year, bene_id, pde_id) %>%
  summarize(antibiotic=unique(antibiotic),
            tier=min(tier)) %>%
  ungroup() %>%
  left_join(bene, by=c('year', 'bene_id'))

pde_approp_table = pde_approp %>%
  count(year, antibiotic, tier) %T>%
  output_table('pde_approp_table')

# inappropriateness (logistic) regression
inapprop_trend = pde_approp %>%
  # outcome is inappropriate
  mutate(y=tier==3) %>%
  glm_f(frmla, family='binomial') %>%
  sandwich_tidy %T>%
  output_table('pde_inapprop_trend')

# what are trends in prescribing practice?
dx_to_pde = read_years(2011:2014, '../../data/dx_to_pde_%i.tsv') %>%
  rename(bene_id=BENE_ID, dx_cat=diagnosis_category) %>%
  # lump abx outside of top 10 into "other"
  mutate(antibiotic=fct_other(factor(antibiotic), keep=top_abx)) %>%
  left_join(bene, by=c('year', 'bene_id'))

# count how many of each encounter type there are in a year
dx_denom = dx_to_pde %>%
  select(year, encounter_id, dx_cat) %>%
  distinct() %>%
  count(year, dx_cat) %>% ungroup() %>%
  rename(n_dx=n)

dx_to_pde_table = dx_to_pde %>%
  count(year, dx_cat, antibiotic) %>% ungroup() %>%
  rename(n_dx_with_abx=n) %>%
  left_join(dx_denom, by=c('year', 'dx_cat')) %>%
  mutate(f_dx_with_abx=n_dx_with_abx/n_dx) %T>%
  output_table('dx_to_pde_table')

encounter = dx_to_pde %>%
  select(year, bene_id, encounter_id, dx_cat) %>%
  distinct()

dx_to_pde_trend_f = function (a, dx) {
  dx_to_pde %>%
    filter(dx_cat==dx) %>%
    group_by(year, encounter_id) %>%
    summarize(bene_id=unique(bene_id), age=unique(age), sex=unique(sex), race=unique(race), dual=unique(dual), n_cc=unique(n_cc), region=unique(region),
              y=a %in% antibiotic) %>%
    ungroup() %>%
    glm_f(frmla, family='binomial') %>%
    sandwich_tidy
}

dx_to_pde_trends = crossing(antibiotic=top_abx, dx_cat=dxs) %>%
  group_by_all() %>%
  do(dx_to_pde_trend_f(.$antibiotic, .$dx_cat)) %>%
  ungroup() %T>%
  output_table('dx_to_pde_trends')

# pde appropriateness at the margins
pde_approp_margin = pde_approp %>%
  mutate(app=tier!=3) %>%
  group_by(year, bene_id, antibiotic, app) %>%
  summarize(n_abx_app_claims=n()) %>%
  group_by(year, bene_id, antibiotic) %>%
  mutate(n_abx_claims=sum(n_abx_app_claims)) %>%
  group_by(year, antibiotic, n_abx_claims, app) %>%
  summarize(n_abx_app_claims=sum(n_abx_app_claims)) %>%
  group_by(year, antibiotic, n_abx_claims) %>%
  mutate(f_app=n_abx_app_claims/sum(n_abx_app_claims)) %>%
  ungroup() %T>%
  output_table('pde_approp_margin')
