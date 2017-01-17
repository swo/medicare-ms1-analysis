#!/usr/bin/env Rscript

bounds_group = function(x, bounds_pairs) {
  out = data_frame(x=x, y='NA')
  for (i in seq_along(bounds_pairs)) {
    bounds = bounds_pairs[[i]]
    if (length(bounds) != 2) {
      print(bounds)
      sprintf("bad bounds: the %i-th bounds has length %i", i, length(bounds)) %>% stop
    }

    lower = bounds[1]
    upper = bounds[2]
    if (lower == -Inf) {
      name = sprintf("(%i) < %i", i, upper)
    } else if(upper == Inf) {
      name = sprintf("(%i) > %i", i, lower)
    } else {
      name = sprintf("(%i) %i -- %i", i, lower, upper)
    }
    out %<>% mutate(y=if_else(between(x, lower, upper), name, y))
  }
  factor(out$y)
}

summarize_by = function(dat, bene, by) {
  dat %>%
    group_by_(by) %>%
    summarize(n_rx=n()) %>%
    left_join(bene %>% group_by_(by) %>% summarize(n_ppl=n()), by=by) %>%
    mutate(rx_per_1k_ppl=n_rx*1000/n_ppl)
}


regions = read_tsv('../../db/census-regions/census-regions.tsv') %>%
  select(state, region)

analyze = function(year) {
  bene_fn = sprintf('../bene-%s.tsv', year)
  pde_fn = sprintf('../abx-pde-%s.tsv', year)

  bene = read_tsv(bene_fn) %>%
    # keep only the beneficiaries that don't have bogus states
    filter(!is.na(state)) %>%
    # merge in the census region information
    left_join(regions, by='state')

  # read in PDE data
  dat = read_tsv(pde_fn) %>%
    # merge the abx class information
    left_join(read_tsv('../abx-classes.tsv'), by='abx') %>%
    # merge the beneficiary information, keeping only claims mapping to beneficiaries
    # for which we have data (i.e., inner join)
    inner_join(bene, by='bene')

  # total prescriptions and prescriptions per 1000 people
  # broken down by antibiotic category
  n_beneficiaries = dim(bene)[1]
  rx_by_class = dat %>%
    group_by(abx_class) %>%
    summarize(n_rx=n()) %>%
    mutate(rx_per_1k_ppl=n_rx*1000/n_beneficiaries)
  write_tsv(rx_by_class, sprintf('rx_by_class-%s.tsv', year))

  # but with a separate list of the top 5 agents
  top_abx_rx = dat %>%
    group_by(abx) %>%
    summarize(n_rx=n()) %>%
    arrange(desc(n_rx)) %>%
    head(5) %>%
    mutate(rx_per_1k_ppl=n_rx*1000/n_beneficiaries)
  write_tsv(top_abx_rx, sprintf('top_abx_rx-%s.tsv', year))

  # from here on, the denominators are taken from the beneficiary data grouped
  # in the same way as the consumption data, so we can use the summarize_by function

  # then broken down by sex
  rx_by_sex = summarize_by(dat, bene, 'sex')
  write_tsv(rx_by_sex, sprintf('rx_by_sex-%s.tsv', year))

  # then by age group
  # (the input to summarize_by is modified, since we want to group by age *group*,
  # not just by raw age)
  age_group_f = function(x) bounds_group(x, list(c(0, 65), c(66, 75), c(76, 85), c(86, Inf)))
  rx_by_age = summarize_by(dat %>% mutate(age_group=age_group_f(age)), bene %>% mutate(age_group=age_group_f(age)), 'age_group')
  write_tsv(rx_by_age, sprintf('rx_by_age-%s.tsv', year))

  # then by census region
  rx_by_region = summarize_by(dat, bene, 'region') %>%
    arrange(desc(rx_per_1k_ppl))
  write_tsv(rx_by_region, sprintf('rx_by_region-%s.tsv', year))

  # prescriptions per 1000 people by state
  # (broken into sextiles to make a colored map)
  rx_by_state = summarize_by(dat, bene, 'state') %>%
    mutate(sextile=ntile(rx_per_1k_ppl, 6))
  write_tsv(rx_by_state, sprintf('rx_by_state-%s.tsv', year))
}

analyze('2011')
analyze('2012')
analyze('2013')
analyze('2014')

# logistic regression on...
# dependent: county in top quartile for prescription rate
# independent:
#  - prescribers per capita above median,
#  - proportion of people obese above median,
#  - proportion of people <2 years old above median,
#  - proportion of people black above median,
#  - per capita income in top 2/3s,
#  - proportion with college education in top 2/3s,
#  - proportion female in top 3/4s
#
# Univariate regressions of each independent against the dependent
# Multivariate regression somehow using adjusted ratios with confounders
