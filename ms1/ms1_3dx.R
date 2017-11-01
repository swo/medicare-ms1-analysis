#!/usr/bin/env Rscript

# load data
# NB: I put DC into the South

dxs = c('gi_t3', 'uti_t3', 'not_infectious', 'uri', 'other_resp', 'gi', 'asthma',
        'ssti_t3', 'uti', 'ssti', 'bronchitis', 'misc_t3', 'om', 'pharyngitis',
        'pneumonia', 'misc_bacterial', 'sinusitis', 'flu', 'om_t3', 'acne', 'viral_flu')

top_abx = c('azithromycin', 'ciprofloxacin', 'amoxicillin', 'cephalexin',
            'trimethoprim/sulfamethoxazole', 'levofloxacin', 'amoxicillin/clavulanate',
            'doxycycline', 'nitrofurantoin', 'clindamycin')

output_table = function(x, base) write_tsv(x, sprintf('tables/%s.tsv', base))

# bene
bene = readRSD('data/bene.rsd')

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

# swo:
# need to be able to count E&M here
dx_from_pde = read_tsv('raw_data/dx_from_pde.tsv')

dx_from_pde %>%
  group_by(year, has_encounter, has_em_encounter) %>%
  summarize(n_pde=sum(n_pde)) %>% ungroup() %>%
  output_table('pde_with_dx_upstream')

# which dx's contribute to each abx?
dx_from_pde %>%
  group_by(year, antibiotic, dx_cat, has_em_encounter) %>%
  summarize(n_pde=sum(n_pde),
            n_pde_with_em_encounter=sum(n_pde[has_em_encounter])) %>%
  ungroup() %>%
  output_table('dx_from_pde')

# appropriateness
fd = read_tsv('db/fd_categories.tsv') %>%
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
dx_to_pde = read_years(2011:2014, 'data/dx_to_pde_%i.tsv') %>%
  left_join(bene, by=c('year', 'bene_id'))

# count how many of each encounter type there are in a year
dx_denom = dx_to_pde %>%
  group_by(year, dx_cat) %>%
  summarize(n_encounters=sum(n_encounters)) %>%
  ungroup()

dx_to_pde_table = dx_to_pde %>%
  group_by(year, dx_cat) %>%
  summarize_at(vars(top_abx), sum) %>%
  ungroup() %>%
  gather('antibiotic', 'n_encounters_with_abx', top_abx) %>%
  left_join(dx_denom, by=c('year', 'dx_cat')) %>%
  mutate(f_encounters_with_abx=n_encounters_with_abx/n_encounters) %T>%
  output_table('dx_to_pde_table')

dxs = c('gi_t3', 'uti_t3', 'uri', 'other_resp', 'gi', 'asthma',
        'ssti_t3', 'uti', 'ssti', 'bronchitis', 'misc_t3', 'om', 'pharyngitis',
        'pneumonia', 'misc_bacterial', 'sinusitis', 'flu', 'om_t3', 'acne', 'viral_flu')

dx_to_pde_trend_f = function (a, dx) {
  dx_to_pde %>%
    filter(dx_cat==dx) %>%
    rename(y=!!a) %>%
    mutate(y=y/n_encounters) %>%
    glm(frmla, data=., family='binomial', weights=n_encounters) %>%
    sandwich_tidy()
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
