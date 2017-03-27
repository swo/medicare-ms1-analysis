bene = read_tsv('../bene_2011.tsv')
pde = read_tsv('../pde_2011.tsv') %>%
  rename(drug=antibiotic) %>%
  group_by(bene_id, drug) %>%
  summarize(n_claims=n()) %>%
  ungroup()

drugs = unique(pde$drug)

# add "all" antibiotics (the sum of usage of all drugs)
#pde %<>%
#  bind_rows(pde %>% group_by(bene_id) %>% summarize(n_claims=sum(n_claims)) %>% mutate(antibiotic='all'))

n_bene = nrow(bene)

user_overlap = function(pde, drug1_name, drug2_name) {
  full_join(pde %>% filter(drug==drug1_name) %>% rename(drug1=drug),
            pde %>% filter(drug==drug2_name) %>% rename(drug2=drug), by='bene_id') %>%
    replace_na(list(drug1=0, drug2=0)) %>%
    mutate(case=case_when(.$drug1 > 0 & .$drug2 > 0 ~ 'took_both',
                          .$drug1 > 0 ~ 'took_drug1_only',
                          .$drug2 > 0 ~ 'took_drug2_only',
                          TRUE ~ 'took_neither')) %>%
    count(case) %>%
    spread(case, n) %>%
    mutate(drug1=drug1_name, drug2=drug2_name)
}

drug_overlap = crossing(drug1=drugs, drug2=drugs) %>%
  filter(drug1 > drug2) %>%
  rowwise() %>% do(user_overlap(pde, .$drug1, .$drug2)) %>%
  mutate(took_neither=n_bene-(took_both + took_drug1_only + took_drug2_only))

pad_to = function(x, n, with=0) {
  if (length(x) > n) stop('input longer than pad to length')
  c(x, rep(0, n - length(x)))
}

top_bracket_limit = function(n, fraction) {
  n[cume_dist(n) > (1 - fraction)] %>% min
}

drug_upper_brackets = pde %>%
  inner_join(pde %>% group_by(bene_id) %>% summarize(n_claims=sum(n_claims)) %>% mutate(antibiotic='all'), by='bene_id') %>%
  group_by(antibiotic) %>%
  summarize(limit=top_bracket_limit(pad_to(n_claims, n_bene), 0.001)) %>%
  filter(limit > 0)  