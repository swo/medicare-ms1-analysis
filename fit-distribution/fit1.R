library(texmexseq)

fit_table = function(n) {
  f = poilogMLE(n, mean(n), sd(n))
  table(n) %>%
    as_data_frame %>%
    mutate(n=as.integer(n)) %>%
    rename(n_claims=n, n_bene=n.1) %>%
    mutate(empirical_d=n_bene/sum(n_bene),
           empirical_cdf=cumsum(empirical_d),
           poilog_d=dpoilog(n_claims, f$par['mu'], f$par['sig']),
           poilog_cdf=cumsum(poilog_d))
}

read_tsv('../pde_2011.tsv') %>%
  rename(drug=antibiotic) %>%
  group_by(bene_id) %>%
  summarize(n_claims=n()) %>%
  count(n_claims) %>%
  write_tsv('n_claims.tsv')