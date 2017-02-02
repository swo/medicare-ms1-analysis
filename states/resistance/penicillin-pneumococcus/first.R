#!/usr/bin/env Rscript

library(broom)

# Do the analysis like before, but just for this one bug/drug combo, and not
# asking if the antibiograms are coherent or anything

bug = 'Streptococcus pneumoniae (non-meningitis)'
consumption_drug=c('amoxicillin', 'penicillin')
resistance_drug=c('ampicillin', 'penicillin')

bug_drug = expand.grid(bug=bug,
                       consumption_drug=consumption_drug,
                       resistance_drug=resistance_drug,
                       stringsAsFactors=FALSE)

# do the correlations

abg = read_tsv('../../../../../antibiogram/data/abg.tsv') %>%
  mutate(drug=tolower(drug)) %>%
  rename(resistance_drug=drug)

known_states = abg %$% state %>% unique %>% sort

state_bene = read_tsv('../../bene-consumption/state-bene-2011.tsv') %>%
  filter(state %in% known_states)
state_cons = read_tsv('../../bene-consumption/state-consumption-2011.tsv') %>%
  filter(state %in% known_states) %>%
  rename(consumption_drug=drug)

sus = bug_drug %>%
  left_join(abg, by=c('bug', 'resistance_drug')) %>%
  group_by(bug, consumption_drug, resistance_drug, state) %>%
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
      sprintf("%s %s %s %s\n", resp, expl, head(df$resistance_drug, 1), head(df$consumption_drug, 1), head(df$bug, 1)) %>% cat
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
  left_join(bug_drug, by=c('bug', 'consumption_drug', 'resistance_drug')) %>%
  left_join(state_cons, by=c('state', 'consumption_drug')) %>%
  group_by(bug, consumption_drug, resistance_drug)

write_tsv(dat, 'plotting-dat.tsv')

res = dat %>%
  do(models(.))

write_tsv(sus, 'susceptibily.tsv')
write_tsv(res, 'model-results.tsv')

res %>%
  filter(p.value < 0.05) %>%
  select(bug, consumption_drug, resistance_drug, explanatory, response, correlation=estimate, p.value) %>%
  mutate(correlation=signif(correlation, digits=3), p.value=signif(p.value, digits=3)) %>%
  write_tsv('model-results-p005.tsv')
