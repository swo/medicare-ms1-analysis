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

read_years = function(years, template) {
  lapply(years, function(y) {
    read_tsv(sprintf(template, y)) %>%
      mutate(year=y)
  }) %>%
    bind_rows()
}

glm_f = function(df, frmla, ...) eval(substitute(function(df2) glm(frmla, data=df2, ...)))(df)
frmla = y ~ year + age + n_cc + sex + race + dual + region

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
  replace_na(list(tier=4)) %>%
  group_by(year, pde_id) %>%
  summarize(antibiotic=unique(antibiotic), tier=min(tier)) %>%
  ungroup() %>%
  left_join(bene, by=c('year', 'bene_id'))

pde_approp_table = pde_approp %>%
  count(year, antibiotic, tier) %T>%
  output_table('pde_approp_table')

# logistic regression
# two ways to do this: either keep the "not infectious" diagnoses in the
# denominator, or leave them out (i.e., filter tier != 4)
# swo: i think this formula is repeated
frmla = y ~ year + age + n_cc + sex + race + dual + region
approp_trends_f = function(df) {
  df %>%
    mutate(y=tier <= 2) %>%
    (function(df) {
      bind_rows(
        glm_f(df, frmla, family='binomial') %>% tidy %>% mutate(antibiotic='overall'),
        df %>%
          filter(antibiotic %in% top_abx) %>%
          group_by(antibiotic) %>%
          do(tidy(glm_f(., frmla, family='binomial'))) %>%
          ungroup()
      )
    })
}

pde_approp_trends = pde_approp %>%
  approp_trends_f() %T>%
  output_table('pde_approp_trends')

# what are trends in prescribing practice?
# swo: test this part
dx_to_pde = read_years(2011:2014, '../../data/dx_to_pde_%i.tsv') %>%
  # hack to fix the weird TMP/SMX values
  mutate(antibiotic=if_else(antibiotic=='trimethoprim/sulfametho',
                            'trimethoprim/sulfamethoxazole',
                            antibiotic)) %>%
  rename(bene_id=BENE_ID, dx_cat=diagnosis_category) %>%
  replace_na(list(antibiotic='no_abx')) %>%
  # drop no infection and no abx, since we won't look at those
  filter(!(dx_cat=='not_infectious' & antibiotic=='no_abx')) %>%
  # lump abx outside of top 10 into "other"
  mutate(antibiotic=fct_other(factor(antibiotic), keep=c(top_abx, 'no_abx'))) %>%
  count(year, bene_id, dx_cat, antibiotic) %>%
  left_join(bene, by=c('year', 'bene_id'))

dx_to_pde_table = dx_to_pde %>%
  group_by(year, dx_cat, antibiotic) %>%
  summarize(n=sum(n)) %>%
  mutate(f=n/sum(n)) %>%
  ungroup() %T>%
  output_table('dx_to_pde_table')

frmla = y ~ year + age + n_cc + sex + race + dual + region
dx_rx_trends = crossing(antibiotic=top_abx, dx_cat=dxs) %>%
  group_by(antibiotic, dx_cat) %>%
  do((function(a, dx) {
    dx_to_pde %>%
      filter(dx_cat==dx) %>%
      mutate(y=antibiotic==a) %>%
      glm_f(frmla, family='binomial', weights=n) %>%
      tidy
  })(.$antibiotic, .$dx_cat)) %>%
  ungroup() %T>%
  output_table('dx_to_pde_trends')

# pde appropriateness at the margins
pde_approp_margin = pde_approp %>%
  # filter out refills if desired!
  mutate(antibiotic=if_else(str_detect(antibiotic, '^trimethoprim/sulfa'), 'tmp_smx', antibiotic)) %>%
  filter(antibiotic %in% c('tmp_smx', top_abx)) %>%
  mutate(app=tier <= 2) %>%
  group_by(year, bene_id, antibiotic, app) %>%
  summarize(n_abx_app_claims=n()) %>%
  group_by(year, bene_id, antibiotic) %>%
  mutate(n_abx_claims=sum(n_abx_app_claims)) %>%
  group_by(year, antibiotic, n_abx_claims, app) %>%
  summarize(n_abx_app_claims=sum(n_abx_app_claims)) %>%
  group_by(year, antibiotic, n_abx_claims) %>%
  mutate(f_app=n_abx_app_claims/sum(n_abx_app_claims)) %>%
  ungroup()
