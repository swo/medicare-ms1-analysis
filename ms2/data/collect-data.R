# inequality function
import::from(ineq, Gini)
nonzero = function(x) x[x > 0]
fraction_nonzero = function(x) length(nonzero(x)) / length(x)

pad0 = function(x, to) c(x, rep(0, to - length(x)))

lm_residuals = function(y, x) residuals(lm(y ~ x))

inequalities = function(x) {
  data_frame(total=sum(x),
             mean=mean(x),
             fnz=fraction_nonzero(x),
             nzgini=Gini(nonzero(x)))
}

consumption_groups = read_tsv('../consumption_groups.tsv')

regions = read_tsv('../../../db/census-regions/census-regions.tsv') %>%
  select(state, state_abbreviation, region)

hrr = read_tsv('../../../db/hrr/hrr.tsv') %>%
  filter(year==2011) %>%
  select(zipcode, hrr)

bene = read_tsv('../../bene_2011.tsv') %>%
  filter(between(age, 66, 100), hmo_months==0) %>%
  left_join(regions, by='state') %>%
  left_join(hrr, by='zipcode')

pde = read_tsv('../../pde_2011.tsv') %>%
  rename(drug=antibiotic) %>%
  left_join(consumption_groups, by='drug') %>%
  group_by(bene_id, drug_group) %>%
  summarize(n_claims=n()) %>%
  ungroup

# compute total claims separately from individual drug groups
total_pde = pde %>%
  group_by(bene_id) %>%
  summarize(n_claims=sum(n_claims)) %>%
  right_join(bene, by='bene_id') %>%
  replace_na(list(n_claims=0)) %>%
  group_by(state) %>%
  do(inequalities(.$n_claims)) %>%
  ungroup() %>%
  mutate(mean_fnz_residual=lm_residuals(mean, fnz))

write_tsv(total_pde, 'all_drug_ineq.tsv')

state_denom = bene %>% count(state) %>% rename(n_bene=n)
state_inequality = pde %>%
  right_join(select(bene, bene_id, state), by='bene_id') %>%
  filter(!is.na(drug_group)) %>%
  left_join(state_denom, by='state') %>%
  group_by(state, drug_group) %>%
  do(inequalities(pad0(.$n_claims, unique(.$n_bene)))) %>%
  ungroup() %>%
  mutate(mean_fnz_residual=lm_residuals(mean, fnz))

write_tsv(state_inequality, 'drug_ineq.tsv')

# repeat this for the HRRs
hrr_total_pde = pde %>%
  group_by(bene_id) %>%
  summarize(n_claims=sum(n_claims)) %>%
  right_join(bene, by='bene_id') %>%
  replace_na(list(n_claims=0)) %>%
  group_by(hrr) %>%
  do(inequalities(.$n_claims)) %>%
  ungroup() %>%
  mutate(mean_fnz_residual=lm_residuals(mean, fnz))

write_tsv(hrr_total_pde, 'hrr_all_drug_ineq.tsv')

hrr_denom = bene %>% count(hrr) %>% rename(n_bene=n)
hrr_inequality = pde %>%
  right_join(select(bene, bene_id, hrr), by='bene_id') %>%
  filter(!is.na(drug_group)) %>%
  left_join(hrr_denom, by='hrr') %>%
  group_by(hrr, drug_group) %>%
  do(inequalities(pad0(.$n_claims, unique(.$n_bene)))) %>%
  ungroup() %>%
  mutate(mean_fnz_residual=lm_residuals(mean, fnz))

write_tsv(hrr_inequality, 'hrr_drug_ineq.tsv')
