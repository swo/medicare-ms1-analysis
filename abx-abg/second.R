# like first, but at the state level

library(purrr)
library(lazyeval)

drugs = c('ciprofloxacin', 'levofloxacin', 'trimethoprim-sulfamethoxazole', 'ampicillin', 'ampicillin-sulbactam')
bugs = c('Escherichia coli')

abx = read_tsv('../bene-abx-2011.tsv') %>%
  # filter out the people who are not old
  filter(age >= 65)

abg = read_tsv('../../../antibiogram/data/abg.tsv', col_types=cols(zipcode=col_character())) %>%
  mutate(drug=tolower(drug))

# compute susceptibility information for each bug/drug combo
sus = abg %>%
  filter(bug %in% bugs, drug %in% drugs) %>%
  group_by(bug, drug, state) %>%
  summarize(n_abg=n(),
            n_isolates=sum(n_isolates, na.rm=TRUE),
            n_years=length(unique(year)),
            median_sus=median(percent_susceptible, na.rm=TRUE),
            mean_sus=mean(percent_susceptible, na.rm=TRUE),
            sd_sus=sd(percent_susceptible, na.rm=TRUE),
            min_sus=min(percent_susceptible, na.rm=TRUE),
            max_sus=max(percent_susceptible, na.rm=TRUE))

write_tsv(sus, 'tmp-sus')

known_states = sus %$% state %>% unique %>% sort

# compute summary stats for the state
state_summary = abx %>%
  filter(state %in% known_states) %>%
  group_by(state) %>%
  summarize(n_bene=n(),
            n_abx_claims=sum(n_claims),
            n_abx_days=sum(days),
            n_zero_abx=sum(days == 0),
            n_pdes=sum(all_pde),
            mean_comorbidity=mean(comorbidity),
            mean_age=mean(age),
            n_white=sum(race=='white'),
            n_male=sum(sex=='male'))

write_tsv(state_summary, 'tmp-hrr')

nz = function(x) x[x > 0]

# compute consumption for each bug/drug combo
drug_vars = map(drugs, as.name)
con = abx %>%
  filter(state %in% known_states) %>%
  select_(.dots=c('bene', 'state', drug_vars)) %>%
  gather_('drug', 'days', drugs) %>%
  group_by(state, drug) %>%
  summarize(n_bene=n(),
            n_days=sum(days),
            n_bene_zero=sum(days==0),
            max_days=max(days),
            mean_nz_days=mean(nz(days)),
            median_nz_days=median(nz(days)),
            quantile_nz_95=quantile(nz(days), probs=c(0.95)),
            quantile_nz_99=quantile(nz(days), probs=c(0.99))) %>%
  mutate(did=1000/365*n_days/n_bene)

write_tsv(con, 'tmp-con')

models = function(df) {
  out = data_frame()
  for (resp in c('median_sus', 'mean_sus', 'min_sus', 'max_sus')) {
    for (expl in c('n_days', 'did', 'n_bene_zero', 'max_days', 'mean_nz_days', 'quantile_nz_95', 'quantile_nz_99')) {
      sprintf("%s %s %s\n", resp, expl, head(df$drug, 1)) %>% cat
      new_rows = df %>%
        mblm(as.formula(paste0(resp, '~', expl)), dataframe=.) %>%
        tidy %>%
        mutate(response=resp, explanatory=expl)
      out = bind_rows(out, new_rows)
    }
  }
  out
}

x = con %>%
  ungroup %>%
  left_join(sus, by=c('drug', 'state')) %>%
  filter(!is.na(n_isolates)) %>%
  arrange(drug)

y = x %>%
  group_by(bug, drug) %>%
  do(models(.))

write_tsv(y, 'tmp-res')
