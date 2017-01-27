library(mblm)
library(broom)

drugs = c('ciprofloxacin', 'levofloxacin', 'trimethoprim-sulfamethoxazole', 'ampicillin', 'ampicillin-sulbactam')
bugs = c('Escherichia coli', 'Coagulase-negative staphylococci')

abg = read_tsv('../../../../antibiogram/data/abg.tsv') %>%
  mutate(drug=tolower(drug))

known_states = abg %$% state %>% unique %>% sort

state_bene = read_tsv('../bene-consumption/state-bene-2011.tsv') %>%
  filter(state %in% known_states)
state_cons = read_tsv('../bene-consumption/state-consumption-2011.tsv') %>%
  filter(state %in% known_states)

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

models = function(df) {
  out = data_frame()
  for (resp in c('median_sus', 'mean_sus', 'min_sus', 'max_sus')) {
    for (expl in c('claims_per_1k_ppl', 'did')) {
      sprintf("%s %s %s\n", resp, expl, head(df$drug, 1), head(df$bug, 1)) %>% cat
      new_rows = cor.test(df[[resp]], df[[expl]], method='spearman') %>%
      #new_rows = df %>%
        #mblm(as.formula(paste0(resp, '~', expl)), dataframe=.) %>%
        tidy %>%
        mutate(response=resp, explanatory=expl)
      out = bind_rows(out, new_rows)
    }
  }
  out
}

x = state_cons %>%
  ungroup %>%
  left_join(sus, by=c('drug', 'state')) %>%
  filter(!is.na(n_isolates)) %>%
  arrange(drug)

y = x %>%
  group_by(bug, drug) %>%
  do(models(.))

write_tsv(y, 'tmp-res')
