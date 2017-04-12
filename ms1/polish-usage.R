regions = read_tsv('../../db/census-regions/census-regions.tsv') %>%
  select(state, region)

antibiotic_class = read_tsv('../abx_class.tsv') %>%
  select(antibiotic, antibiotic_class) %>%
  distinct

analyze = function(year) {
  cc = read_tsv(sprintf('../cc_%i.tsv', year)) %>%
    select(-ALZH)

  n_cc = cc %>% select(AMI:HYPOTH) %>% rowSums

  bene = read_tsv(sprintf('../bene_%i.tsv', year)) %>%
    filter(between(age, 66, 100), hmo_months==0) %>%
    mutate(is_female=sex=='female', is_white=race=='white') %>%
    left_join(regions, by='state') %>%
    left_join(cc, by='bene_id') %>%
    select(-sex, -race, -zipcode, -hmo_months)

  n_cc = bene %>% select(AMI:HYPOTH) %>% rowSums

  bene %<>% mutate(n_cc=n_cc)

  pde = read_tsv(sprintf('../pde_%i.tsv', year)) %>%
    semi_join(bene, by='bene_id') %>%
    mutate(year=year)

  usage = pde %>% count(bene_id) %>% rename(n_claims=n) %>%
    right_join(bene, by='bene_id') %>%
    replace_na(list(n_claims=0)) %>%
    mutate(year=year)

  write_tsv(usage, sprintf('usage_%i.tsv', year))
  write_tsv(pde, sprintf('pde_%i.tsv', year))
}

for (y in 2011:2014) analyze(y)
