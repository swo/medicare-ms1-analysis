#!/usr/bin/env Rscript

# Polish the cost-use data, keeping only beneficiaries that are in the bene data.

polish_cu = function(year) {
  bene_fn = paste0('bene_', year, '.tsv')
  data_fn = paste0('../data/cu_', year, '.tsv')
  output_fn = paste0('cu_', year, '.tsv')

  bene_ids = read_tsv(bene_fn) %$% bene_id

  x = read_tsv(data_fn) %>%
    rename(bene_id=BENE_ID,
          hospital_op_visits=HOP_VISI,
          hospital_er_visits=HOP_ER_V,
          acute_ip_stays=ACUTE_ST,
          snf_stays=SNF_STAY,
          pde=PTD_EVEN) %>%
    filter(bene_id %in% bene_ids)

  # Get the names to replace
  ns = x %>% select(-bene_id) %>% names
  os = rep(0, length(ns))
  rl = setNames(os, ns) %>% as.list

  x %<>% replace_na(rl)

  write_tsv(x, output_fn)
}

for (year in 2011:2014) polish_cu(year)
