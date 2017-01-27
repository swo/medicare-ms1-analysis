#!/usr/bin/env Rscript

coherence = read_tsv('../../../../antibiogram/analysis/geographic-coherence/coherences.tsv') %>%
  rename(resistance_drug=drug) %>%
  mutate(resistance_drug=tolower(resistance_drug)) %>%
  select(bug, resistance_drug, level, n_groups)

gram_negs = c('Escherichia coli', 'Klebsiella pneumoniae', 'Proteus mirabilis',
  'Enterobacter cloacae', 'Klebisella oxytoca')
gram_poss = c('Coagulase-negative staphylococci', 'Methicillin-resistant S. aureus (MRSA)',
  'Methicillin-susceptible S. aureus (MSSA)', 'Streptococcus pneumoniae (non-meningitis)',
  'Staphylococcus aureus')
organisms = c(gram_negs, gram_poss)

cons_res = data_frame(
  consumption_drug=c('trimethoprim-sulfamethoxazole', 'ciprofloxacin', 'levofloxacin', 'ceftriaxone', 'nitrofurantoin', 'amoxicillin', 'penicillin', 'doxycycline', 'azithromycin'),
  resistance_drug=c('trimethoprim-sulfamethoxazole', 'ciprofloxacin', 'levofloxacin', 'ceftriaxone', 'nitrofurantoin', 'ampicillin', 'penicillin', 'tetracycline/doxycycline', 'erythromycin')
)

add_organisms = function(x) {
  data_frame(consumption_drug=x$consumption_drug,
             resistance_drug=x$resistance_drug,
             bug=organisms)
}

# get all combos of cons_res with each bug
# keep only those res-drug/bug combos are geographically coherent
coherent_bug_drug = cons_res %>%
  rowwise %>%
  do(add_organisms(.)) %>%
  ungroup %>%
  left_join(coherence, by=c('resistance_drug', 'bug')) %>%
  mutate(state_coherent=level %in% c('state', 'hrr', 'zipcode'))

# do the correlations

abg = read_tsv('../../../../antibiogram/data/abg.tsv') %>%
  mutate(drug=tolower(drug)) %>%
  rename(resistance_drug=drug)

known_states = abg %$% state %>% unique %>% sort

state_bene = read_tsv('../bene-consumption/state-bene-2011.tsv') %>%
  filter(state %in% known_states)
state_cons = read_tsv('../bene-consumption/state-consumption-2011.tsv') %>%
  filter(state %in% known_states) %>%
  rename(consumption_drug=drug)

sus = coherent_bug_drug %>%
  left_join(abg, by=c('bug', 'resistance_drug')) %>%
  group_by(bug, resistance_drug, state) %>%
  summarize(n_abg=n(),
            n_isolates=sum(n_isolates, na.rm=TRUE),
            n_years=length(unique(year)),
            median_sus=median(percent_susceptible, na.rm=TRUE),
            mean_sus=mean(percent_susceptible, na.rm=TRUE),
            sd_sus=sd(percent_susceptible, na.rm=TRUE),
            min_sus=min(percent_susceptible, na.rm=TRUE),
            max_sus=max(percent_susceptible, na.rm=TRUE))

models = function(df) {
  out = data_frame()
  for (resp in c('median_sus', 'mean_sus', 'min_sus', 'max_sus')) {
    for (expl in c('claims_per_1k_ppl', 'did')) {
      sprintf("%s %s %s\n", resp, expl, head(df$resistance_drug, 1), head(df$bug, 1)) %>% cat
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

dat = sus %>%
  left_join(coherent_bug_drug, by=c('bug', 'resistance_drug')) %>%
  left_join(state_cons, by=c('state', 'consumption_drug')) %>%
  filter(state_coherent) %>%
  group_by(bug, resistance_drug)

write_tsv(dat, 'plotting-dat.tsv')

res = dat %>%
  do(models(.))

write_tsv(sus, 'susceptibily.tsv')
write_tsv(res, 'model-results.tsv')

res %>%
  filter(p.value < 0.05) %>%
  select(bug, resistance_drug, explanatory, response, correlation=estimate, p.value) %>%
  mutate(correlation=signif(correlation, digits=3), p.value=signif(p.value, digits=3)) %>%
  write_tsv('model-results-p005.tsv')
