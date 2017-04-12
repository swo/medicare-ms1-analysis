import::from(lazyeval, interp)
import::from(broom, tidy)
import::from(stringr, str_replace)
import::from(doMC, registerDoMC)
import::from(gplots, heatmap.2)
library(glmnet)

# load db data
regions = read_tsv('../../../db/census-regions/census-regions.tsv') %>%
  select(state, region)
condition_names = read_tsv('../../chronic-conditions/names.txt',
                           col_names=c('condition', 'condition_name'))

# load only data from 2011
usage = read_tsv('../usage_2011.tsv') %>%
  mutate(age_over_67 = age - 67)

conditions = c("AMI", "ALZHDMTA", "ATRIALFB", "CATARACT", "CHRNKIDN", "COPD", "CHF", "DIABETES", "GLAUCOMA", "HIPFRAC", "ISCHMCHT", "DEPRESSN", "OSTEOPRS", "RA_OA", "STRKETIA", "CNCRBRST", "CNCRCLRC", "CNCRPRST", "CNCRLUNG", "CNCRENDM", "ANEMIA", "ASTHMA", "HYPERL", "HYPERP", "HYPERT", "HYPOTH")

# mcc functions
mcc = function(true, pred) {
  ut = true - mean(true)
  up = pred - mean(pred)

  co = mean(ut * up)
  vt = mean(ut ** 2)
  vp = mean(up ** 2)

  if (vt * vp == 0) {
    NA
  } else {
    co / sqrt(vt * vp)
  }
}

mcc_cor = function(x) {
  x = as.matrix(x)
  labels = colnames(x)
  n = length(labels)
  out = matrix(nrow=n, ncol=n, dimnames=list(labels, labels))

  for (i in 1:n) out[i, i] = 1.0

  for (i in 1:(n-1)) {
    ui = x[,i] - mean(x[,i])
    vi = mean(ui ** 2)
    for (j in (i+1):n) {
      uj = x[,j] - mean(x[,j])
      vj = mean(uj ** 2)
      co = mean(ui * uj)
      val = co / sqrt(vi * vj)
      out[i, j] = val
      out[j, i] = val
    }
  }
  out
}

# All-condition model
multivariate_formula = paste(c('n_claims ~ age_over_67 + is_female + ', conditions), collapse=' + ') %>% formula
start = rep(1.0, 3 + length(conditions))

all_condition_model = usage %>%
  filter(age >= 67) %>%
  glm(multivariate_formula, data=., family=quasipoisson(link='identity'), start=start) %>%
  tidy()

write_tsv(all_condition_model, 'table_all_condition_raw.tsv')

# Polished version of that table
ccnames = read_tsv('../../chronic-conditions/names.txt', col_names=c('term', 'condition'))

all_condition_model %>%
  mutate(term=str_replace(term, 'TRUE$', '')) %>%
  left_join(ccnames, by='term') %>%
  mutate(term=if_else(is.na(condition), term, condition)) %>%
  mutate(estimate=round(estimate, digits=3)) %>%
  select(term, estimate) %>%
  write_tsv('table_all_condition.tsv')

# All-condition Poisson lasso model
registerDoMC(cores=3)

lasso_x = usage %>%
  select(age_over_67, is_female, AMI:HYPOTH) %>%
  as.matrix
lasso_y = usage$n_claims

lasso_cv = cv.glmnet(lasso_x, lasso_y, parallel=TRUE, family='poisson')

pdf('supp/lasso_cv.pdf')
plot(lasso_cv)
dev.off()

lasso_df = lasso_cv$glmnet.fit %>%
  tidy %>%
  filter(estimate > 0) %>%
  group_by(step) %>%
  summarize(df=n() - 1)

lasso_cv$glmnet.fit %>%
  tidy %>%
  select(step, lambda, dev.ratio) %>%
  distinct() %>%
  left_join(lasso_df, by='step') %>%
  select(step, df, lambda) %>%
  mutate_at(vars(lambda), function(x) signif(x, digits=2)) %>%
  #filter(df < 12) %>%
  write_tsv('supp/lasso_cv_df.tsv')

lambda_lasso = 0.22

lasso_cv %>%
  coef(s=lambda_lasso) %>%
  tidy %>%
  select(condition=row, coefficient=value) %>%
  left_join(condition_names, by='condition') %>%
  mutate(term=if_else(is.na(condition_name), condition, condition_name)) %>%
  select(term, coefficient) %>%
  arrange(desc(coefficient)) %>%
  write_tsv('supp/lasso_coefficients.tsv')

# Single-condition Poisson models including race and geography for top conditions
top_conditions = c('COPD', 'ASTHMA', 'CNCRLUNG', 'CHF', 'DEPRESSN')

race_geography_model_f = function(condition) {
  form = formula(paste0('n_claims ~ age_over_67 + is_female + is_white + region + ', condition))
  usage %>%
    filter(age >= 67) %>%
    mutate(region=factor(region, levels=c('Northeast', 'West', 'Midwest', 'South'))) %>%
    glm(form, family=quasipoisson(link='identity'), data=.) %>%
    tidy() %>%
    mutate(condition=condition)
}

lapply(top_conditions, race_geography_model_f) %>%
  bind_rows() %>%
  select(condition, term, estimate) %>%
  mutate(term=if_else(term==paste0(condition, 'TRUE'), 'condition', term)) %>%
  spread(term, estimate) %>%
  write_tsv('supp/georace.tsv')

# Correlations between chronic conditions
correlations = usage %>%
  select(AMI:HYPOTH) %>%
  mcc_cor

# Heatmap. Sqrt of MCC used. Female conditions are excluded.
female_conditions = c('CNCRENDM', 'CNCRBRST', 'OSTEOPRS')
pdf('supp/mcc.pdf')
(!(rownames(correlations) %in% female_conditions)) %>%
  { correlations[., .] } %>%
  sqrt %>%
  gplots::heatmap.2(trace='none', symm=TRUE)
dev.off()
