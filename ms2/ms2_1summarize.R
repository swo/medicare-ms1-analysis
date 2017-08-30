#!/usr/bin/env Rscript

library(rlang)

# read in HRR designations for all years
hrr_all = read_tsv('../../db/hrr/hrr.tsv') %>%
  mutate(hrr=as.character(hrr))

consumption_groups = read_tsv('consumption_groups.tsv')

summarize_inequality = function(y) {
  # get HRR designations just for that year
  hrr = hrr_all %>%
    filter(year==y) %>%
    select(zipcode, hrr)

  # set up input/output filenames
  bene_fn = sprintf('../bene_%i.tsv', y)
  pde_fn = sprintf('../pde_%i.tsv', y)

  # get information about bene locations
  bene = read_tsv(bene_fn) %>%
    left_join(hrr, by='zipcode') %>%
    select(bene_id, state, hrr)

  # get PDEs
  pde = read_tsv(pde_fn) %>%
    rename(drug=antibiotic) %>%
    inner_join(consumption_groups, by='drug') %>%
    count(bene_id, drug_group) %>% rename(n_claims=n) %>%
    # add the "overall" category
    (function(df) {
      group_by(df, bene_id) %>%
        summarize(n_claims=sum(n_claims)) %>%
        mutate(drug_group='overall') %>%
        bind_rows(df)
    }) %>%
    left_join(bene, by='bene_id')

  # function to compute the denominators, i.e., the number of benes in each unit
  denom_f = function(bene, unit_type_) {
    bene %>%
      rename(unit=!!sym(unit_type_)) %>%
      count(unit) %>% rename(n_bene=n)
  }

  # function to compute inequality, i.e., the number of beneficiaries taking
  # each number of drug
  ineq_f = function(pde, bene, unit_type_) {
    # get denominator, which we'll use to fill in the non-consumers
    denom = denom_f(bene, unit_type_)

    pde %>%
      # rename, e.g., the "state" column to "unit"
      rename(unit=!!sym(unit_type_)) %>%
      # get, e.g., Massachusetts, amoxicillin, 2 claims, 234 beneficiaries
      count(unit, drug_group, n_claims) %>% rename(n_bene=n) %>%
      # find the zero-consumers
      (function(df) {
        # count un beneficiaries who took >= 1 drug in each unit/drug
        group_by(df, unit, drug_group) %>%
          summarize(n_consumers=sum(n_bene)) %>%
          ungroup() %>%
          # get the denominators
          left_join(denom, by='unit') %>%
          # infer number of nonconsumers from number of consumers and denominator
          mutate(n_nonconsumers = n_bene - n_consumers,
                 n_claims = 0) %>%
          select(unit, drug_group, n_claims, n_bene=n_nonconsumers) %>%
          # add these non-consumer rows to the rest of the data
          bind_rows(df)
      }) %>%
      arrange(unit, drug_group) %>%
      mutate(unit_type=unit_type_)
  }

  # get inequalities for each drug/unit at the state and HRR levels
  ineq = bind_rows(
    ineq_f(pde, bene, 'state'),
    ineq_f(pde, bene, 'hrr')
  ) %>%
    mutate(year=y)

  ineq
}

# get inequalities for each drug, unit, year, and level
lapply(2011:2014, summarize_inequality) %>%
  bind_rows() %>%
  write_tsv('ineq.tsv')
