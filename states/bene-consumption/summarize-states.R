#!/usr/bin/env Rscript

f_state_bene_summary = function(bene) {
  # chronic conditions are missing, they should come here
  bene %>%
    mutate(is_female=if_else(.$sex=='female', 1, 0),
           is_white=if_else(.$race=='white', 1, 0)) %>%
    group_by(state) %>%
    summarize(n_bene=n(),
              fraction_white=sum(is_white)/n_bene,
              fraction_female=sum(is_female)/n_bene,
              fraction_geq_65=sum(age >= 65)/n_bene,
              mean_age=mean(age))
}

f_state_consumption = function(abx, bene, state_bene_summary) {
  abx %>%
    left_join(select(bene, bene, state), by='bene') %>%
    select(state, drug=abx, n_claims, days) %>%
    group_by(state, drug) %>%
    summarize(claims=sum(n_claims), days=sum(days)) %>%
    left_join(select(state_bene_summary, state, n_bene), by='state') %>%
    mutate(claims_per_1k_ppl=claims/n_bene * 1000,
           days_per_1k_ppl=days/n_bene * 1000) %>%
    select(-n_bene)
}

analyze = function(year) {
  bene_fn = sprintf('../../bene-%s.tsv', year)
  abx_fn = sprintf('../../abx-pde-%s.tsv', year)
  state_bene_summary_fn = sprintf('state-bene-%s.tsv', year)
  state_consumption_fn = sprintf('state-consumption-%s.tsv', year)

  bene = read_tsv(bene_fn) %>%
    filter(!is.na(state))

  state_bene_summary = f_state_bene_summary(bene)

  abx = read_tsv(abx_fn)
  state_consumption = f_state_consumption(abx, bene, state_bene_summary)

  write_tsv(state_bene_summary, state_bene_summary_fn)
  write_tsv(state_consumption, state_consumption_fn)
}

analyze('2011')
analyze('2012')
analyze('2013')
analyze('2014')
