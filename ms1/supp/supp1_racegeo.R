import::from(broom, tidy)

# load db data
regions = read_tsv('../../../db/census-regions/census-regions.tsv') %>%
  select(state, region)
condition_names = read_tsv('../../chronic-conditions/names.txt',
                           col_names=c('condition', 'condition_name'))

# load only data from 2011
usage = read_tsv('../usage_2011.tsv') %>%
  mutate(age_over_67 = age - 67)

conditions = c("AMI", "ALZHDMTA", "ATRIALFB", "CATARACT", "CHRNKIDN", "COPD", "CHF", "DIABETES", "GLAUCOMA", "HIPFRAC", "ISCHMCHT", "DEPRESSN", "OSTEOPRS", "RA_OA", "STRKETIA", "CNCRBRST", "CNCRCLRC", "CNCRPRST", "CNCRLUNG", "CNCRENDM", "ANEMIA", "ASTHMA", "HYPERL", "HYPERP", "HYPERT", "HYPOTH")

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

transpose_tibble = function(x) {
  gather(x, 'var_name', 'value', 2:ncol(x)) %>% 
    spread_(names(x)[1], 'value')
}

lapply(top_conditions, race_geography_model_f) %>%
  bind_rows() %>%
  mutate(estimate=round(estimate, digits=3)) %>%
  select(condition, term, estimate) %>%
  mutate(term=if_else(term==paste0(condition, 'TRUE'), 'condition', term)) %>%
  spread(term, estimate) %>%
  write_tsv('table_georace.tsv')
