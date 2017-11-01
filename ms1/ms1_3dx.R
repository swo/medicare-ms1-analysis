#!/usr/bin/env Rscript

output_table = function(x, base) write_tsv(x, sprintf('tables/%s.tsv', base))

# bene
regions = read_tsv('db/census-regions.tsv') %>%
  select(state, region)

bene = read_tsv('raw_data/bene.tsv') %>%
  rename(bene_id=BENE_ID) %>%
  left_join(regions, by='state')

glm_f = function(df, frmla, ...) eval(substitute(function(df2) glm(frmla, data=df2, ...)))(df)
sandwich_tidy = function(m) tidy(lmtest::coeftest(m, vcov=sandwich::vcovHC(m, type='HC3')))

# appropriateness trend
pde_inapprop = read_tsv('raw_data/bene_inapprop.tsv') %>%
  filter(n_firstfill>0) %>%
  mutate(risk=n_inapprop_firstfill/n_firstfill) %>%
  rename(bene_id=BENE_ID) %>%
  left_join(bene, by=c('year', 'bene_id'))

inapprop_risk_2011 = pde_inapprop %>%
  filter(year==2011) %$%
  mean(risk) %T>%
  write('tables/inapprop_risk_2011.txt')

inapprop_trend = pde_inapprop %>%
  glm_f(risk ~ year + age + n_cc + sex + race + dual + region, weights=n_firstfill, family='binomial') %>%
  sandwich_tidy() %T>%
  output_table('inapprop_trend')

# trends in prescribing practice
dx_to_pde = read_tsv('raw_data/dx_to_pde.tsv') %>%
  rename(bene_id=BENE_ID) %>%
  left_join(bene, by=c('year', 'bene_id'))

dx_to_pde_table = dx_to_pde %>%
  select(year, dx_cat, azithromycin:clindamycin) %>%
  mutate(n_encounters=1) %>%
  group_by(year, dx_cat) %>%
  summarize_all(sum) %>%
  ungroup() %>%
  mutate_at(vars(azithromycin:clindamycin), function(x) x / .$n_encounters) %T>%
  output_table('dx_to_pde_table')

dxs = c('gi_t3', 'uti_t3', 'uri', 'other_resp', 'gi',
        'asthma', 'ssti_t3', 'uti', 'ssti', 'bronchitis', 'misc_t3', 'om',
        'pharyngitis', 'pneumonia', 'misc_bacterial', 'sinusitis', 'flu',
        'om_t3', 'acne', 'viral_flu')

abx = c('azithromycin', 'ciprofloxacin', 'amoxicillin', 'cephalexin',
            'tmpsmx', 'levofloxacin', 'amoxclav',
            'doxycycline', 'nitrofurantoin', 'clindamycin')

dx_to_pde_trend_f = function (a, dx) {
  frmla = str_interp("${a} ~ year + age + n_cc + sex + race + dual + region") %>%
    as.formula
  
  now = format(Sys.time(), "%H:%M.%S")
  cat(str_interp("${now} -- abx: ${a} dx: ${dx}\n"))
  
  dx_to_pde %>%
    filter(dx_cat==dx) %>%
    glm_f(frmla, family='binomial') %>%
    sandwich_tidy()
}

dx_to_pde_trends = crossing(abx=abx, dx_cat=dxs) %>%
  group_by_all() %>%
  do(dx_to_pde_trend_f(.$abx, .$dx_cat)) %>%
  ungroup() %T>%
  output_table('dx_to_pde_trends')
