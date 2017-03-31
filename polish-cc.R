#!/usr/bin/env Rscript

# Polish the chronic condition data, relying on the bene data already being done
# swo> Merge this into the other polishing script

polish_cc = function(x) {
  rename(x, bene_id=BENE_ID) %>%
    select(-STRCTFLG) %>%
    mutate_at(vars(AMI:HYPOTH), function(x) x==3)
}

polish_year = function(y) {
  cc_fn = sprintf('../data/cc_%i.tsv', y)
  bene_fn = sprintf('bene_%i.tsv', y)
  cc_out_fn = sprintf('cc_%i.tsv', y)

  bene_ids = read_tsv(bene_fn) %$% bene_id

  cc = read_tsv(cc_fn) %>%
    polish_cc() %>%
    filter(bene_id %in% bene_ids)

  write_tsv(cc, cc_out_fn)
}

for (year in 2011:2014) polish_year(year)
