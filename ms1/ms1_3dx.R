library(forcats)

# load data
# NB: I put DC into the South

output_table = function(x, base) write_tsv(x, sprintf('tables/%s.tsv', base))

# bene
bene = read_tsv('data/bene.tsv') %>%
  mutate(region=factor(region, levels=c('Northeast', 'West', 'Midwest', 'South')),
         sex=factor(sex, levels=c('male', 'female')),
         race=factor(race, levels=c('white', 'black', 'Hispanic', 'other'))) %>%
  mutate(age=age-1)

# trends in diagnoses by beneficiary
read_years = function(years, template) {
  lapply(years, function(y) {
    read_tsv(sprintf(template, y)) %>%
      mutate(year=y)
  }) %>%
    bind_rows()
}

dx_by_bene = read_years(2011:2014, '../../data/dx_by_bene_%i.tsv') %>%
  rename(bene_id=BENE_ID, dx_cat=diagnosis_category)

dxs = unique(dx_by_bene$dx_cat)

glm_f = function(df, frmla, ...) eval(substitute(function(df2) glm(frmla, data=df2, ...)))(df)
frmla = y ~ year + age + n_cc + sex + race + dual + region

dx_trends_models = lapply(dxs, function(dx) {
  dx_by_bene %>%
    filter(dx_cat==dx) %>%
    right_join(bene, by=c('year', 'bene_id')) %>%
    replace_na(list(n_dx=0)) %>%
    rename(y=n_dx) %>%
    glm_f(frmla, family='poisson') %>%
    tidy %>%
    mutate(model=dx)
}) %>%
  bind_rows() %T>%
  output_table('dx_trends')

# which dx's contribute to each abx?
# in 2011
dx_from_pde = read_tsv(sprintf('../../data/dx_from_pde_%i.tsv', 2011)) %>%
  rename(pde_id=PDE_ID, dx_cat=diagnosis_category)

top_abx = c('azithromycin', 'ciprofloxacin', 'amoxicillin', 'cephalexin',
            'trimethoprim/sulfamethoxazole', 'levofloxacin', 'amoxicillin/clavulanate',
            'doxycycline', 'nitrofurantoin', 'clindamycin')

# get a denominator for total number of PDEs
pde_denom = dx_from_pde %>%
  filter(antibiotic %in% top_abx) %>%
  select(pde_id, antibiotic) %>%
  distinct()

dx_from_pde_table = crossing(antibiotic=top_abx, dx_cat=dxs) %>%
  group_by(antibiotic, dx_cat) %>%
  do((function(a, dx) {
    # how what fraction of PDEs for this abx had this dx upstream?
    pde_denom %>%
      filter(antibiotic==a) %>%
      left_join(filter(dx_from_pde, dx_cat==dx), by='pde_id') %>%
      mutate(present=!is.na(dx_cat)) %>%
      count(present)
  })(.$antibiotic, .$dx_cat)) %T>%
  output_table('dx_from_pde')

# what are trends in prescribing practice?
# swo: test this part
dx_to_pde = read_years(2011:2014, '../../data/dx_to_pde_%i.tsv') %>%
  rename(bene_id=BENE_ID, dx_cat=diagnosis_category) %>%
  mutate(dummy=TRUE) %>% spread(antibiotic, dummy) %>%
  replace_na(as.list(rep(0, length(top_abx))), top_abx) %>%
  left_join(bene, by=c('year', 'bene_id'))

dx_rx_table = crossing(antibiotic=top_abx, dx_cat=dxs) %>%
  group_by(antibiotic, dx_cat) %>%
  do((function(a, dx) {
    dx_to_pde %>%
      filter(dx_cat==dx) %>%
      count(!!sym(a))
  })(.$antibiotic, .$dx_cat)) %T>%
  output_table('dx_rx_table')

dx_rx_trends = crossing(antibiotic=top_abx, dx_cat=dxs) %>%
  group_by(antibiotic, dx_cat) %>%
  do((function(a, dx) {
    frmla = as.formula(str_interp("${a} ~ year + age + n_cc + sex + race + dual + region"))
    dx_to_pde %>%
      filter(dx_cat==dx) %>%
      glm_f(frmla, family='binomial') %>%
      tidy
  }) %T>%
  output_table('dx_rx_trends')
