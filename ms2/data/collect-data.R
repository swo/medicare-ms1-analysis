# inequality function
import::from(ineq, Gini)
nonzero = function(x) x[x > 0]
fraction_nonzero = function(x) length(nonzero(x)) / length(x)

pad0 = function(x, to) c(x, rep(0, to - length(x)))

inequalities = function(x) {
  data_frame(total=sum(x),
             mean=mean(x),
             fnz=fraction_nonzero(x),
             nzgini=Gini(nonzero(x)))
}

consumption_groups = read_tsv('../consumption_groups.tsv')

regions = read_tsv('../../../db/census-regions/census-regions.tsv') %>%
  select(state, state_abbreviation, region)

bene = read_tsv('../../bene_2011.tsv') %>%
  filter(between(age, 66, 100), hmo_months==0) %>%
  left_join(regions, by='state')


# compute total claims separately from individual drug groups
total_pde = pde %>%
  group_by(bene_id) %>%
  summarize(n_claims=sum(n_claims)) %>%
  right_join(bene, by='bene_id') %>%
  replace_na(list(n_claims=0)) %>%
  group_by(state) %>%
  do(inequalities(.$n_claims))

write_tsv(total_pde, 'total_pde.tsv')

state_denom = bene %>% count(state) %>% rename(n_bene=n)
state_inequality = pde %>%
  right_join(select(bene, bene_id, state), by='bene_id') %>%
  filter(!is.na(drug_group)) %>%
  left_join(state_denom, by='state') %>%
  group_by(state, drug_group) %>%
  do(inequalities(pad0(.$n_claims, unique(.$n_bene))))

write_tsv(state_inequality, 'group_inequality.tsv')
