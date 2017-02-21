bbx = read_tsv('bene-abx-2011.tsv')

bcor = bbx %>%
  select(-bene, -n_claims, -total_abx, -all_pde, -age, -zip, -plan_coverage_months, -county, -state, -race, -comorbidity) %>%
  cor(method='spearman') %>%
  as.data.frame %>%
  mutate(abx1=names(.)) %>%
  gather(abx2, spearman, amikacin:vancomycin) %>%
  filter(abx1 < abx2)

sample_size = 10
bcor.sample = sample_n(bcor, sample_size)
pvalues = rep(NA, length.out=sample_size)

for (i in 1:sample_size) {
  pvalues[i] = cor.test(bbx[[bcor$abx1[i]]], bbx[[bcor$abx2[i]]], method='spearman', exact=FALSE)$p.value
}

bcor.sample$pvalues = pvalues