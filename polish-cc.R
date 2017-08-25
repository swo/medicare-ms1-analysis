#!/usr/bin/env Rscript

# Polish the chronic condition data, relying on the bene data already being done
# swo> Merge this into the other polishing script

polish_cc = function(x) {
  rename(x, bene_id=BENE_ID) %>%
    mutate_at(vars(AMI:HYPOTH), function(x) x==3)
}

polish_year = function(y) {
  cc_fn = sprintf('../data/cc_%i.tsv', y)
  bene_fn = sprintf('bene_%i.tsv', y)
  cc_out_fn = sprintf('cc_%i.feather', y)

  bene = read_tsv(bene_fn)

  read_tsv(cc_fn) %>%
    polish_cc() %>%
    semi_join(bene, by='bene_id') %>%
    write_feather(cc_out_fn)
}

library(parallel)
cluster = makeCluster(3, type='FORK')

parLapply(cluster, 2011:2014, polish_year)

stopCluster(cluster)
