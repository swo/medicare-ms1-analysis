#!/usr/bin/env Rscript

# Summarize the information about abx usage for each beneficiary:
# - Use the abx PDE information
# - Give just # of events, total days, and total cost

summarize_abx_pde = function(in_fn, out_fn) {
  read_tsv(in_fn) %>%
    rename(bene=BENE_ID) %>%
    group_by(bene) %>%
    summarize(total_abx_pde=n(), total_days=sum(DAYSSPLY), total_cost=sum(TOTALCST)) %>%
    select(bene, total_abx_pde, total_days, total_cost) %>%
    write_tsv(out_fn)
}

summarize_abx_pde('../data/abx_pde_2011.tsv', 'summary-abx-pde-2011.tsv')