#!/usr/bin/env Rscript

# utility functions function
import::from(ineq, ineq)
nonzero = function(x) x[x > 0]
fraction_nonzero = function(x) length(nonzero(x)) / length(x)
pad0 = function(x, to) c(x, rep(0, to - length(x)))
filter_eq_ = function(df, col, val) filter_(df, str_interp("${col}=='${val}'"))

inequalities = function(pde, drug_group, unit_type, unit, n_bene) {
  if (drug_group=='overall') {
    filter_drug_f = function(df) df
  } else {
    filter_drug_f = function(df) filter_eq_(df, 'drug_group', drug_group)
  }

  x = pde %>%
    filter_eq_(unit_type, unit) %>%
    filter_drug_f() %$%
    pad0(n_claims, n_bene)

  nb_par = fitdistrplus::mledist(x, 'nbinom')$estimate
  nb_size = nb_par['size']
  nb_prob = nb_size / (nb_size + nb_par['mu'])

  data_frame(drug_group=drug_group,
             unit_type=unit_type,
             unit=unit,
             n_bene=n_bene,
             n=length(x),
             total=sum(x),
             mean=mean(x),
             fnz=fraction_nonzero(x),
             gini=ineq(x, type='Gini'),
             nzgini=ineq(nonzero(x), type='Gini'),
             nb_size=nb_size,
             nb_prob=nb_prob)
}

# auxiliary data
consumption_groups = read_tsv('consumption_groups.tsv')

regions = read_tsv('../../db/census-regions/census-regions.tsv') %>%
  select(state, state_abbreviation, region) %>%
  filter(state != 'District of Columbia')
stopifnot('DC' %in% regions$state)

# read in HRR designations for all years
hrr_all = read_tsv('../../db/hrr/hrr.tsv') %>%
  mutate(hrr=as.character(hrr))

summarize_inequality = function(y) {
  # get HRR designations just for that year
  hrr = hrr_all %>%
    filter(year==y) %>%
    select(zipcode, hrr)

  # set up input/output filenames
  bene_fn = sprintf('../bene_%i.tsv', y)
  pde_fn = sprintf('../pde_%i.tsv', y)

  bene = read_tsv(bene_fn) %>%
    left_join(regions, by='state') %>%
    left_join(hrr, by='zipcode') %>%
    select(bene_id, state, zipcode, region, hrr)

  pde = read_tsv(pde_fn) %>%
    rename(drug=antibiotic) %>%
    left_join(consumption_groups, by='drug') %>%
    group_by(bene_id, drug_group) %>%
    summarize(n_claims=n()) %>%
    ungroup() %>%
    left_join(bene, by='bene_id')

  # compute number of benes in each state, HRR
  denom = bind_rows(
    bene %>% rename(unit=state) %>% count(unit) %>% mutate(unit_type='state'),
    bene %>% rename(unit=hrr) %>% count(unit) %>% mutate(unit_type='hrr')
  ) %>%
    rename(n_bene=n)

  # compute and return inequalities
  crossing(drug_group=c('overall', consumption_groups$drug_group), unit=denom$unit) %>%
    left_join(denom, by='unit') %>%
    rowwise() %>%
    do(inequalities(pde, .$drug_group, .$unit_type, .$unit, .$n_bene)) %>%
    ungroup() %>%
    mutate(year=y)
}

lapply(2011:2014, summarize_inequality) %>%
  bind_rows() %>%
  write_tsv('ineq.tsv')
