#!/usr/bin/env Rscript

coherence = read_tsv('../../../../antibiogram/analysis/geographic-coherence/coherences.tsv') %>%
  rename(resistance_drug=drug) %>%
  mutate(resistance_drug=tolower(resistance_drug)) %>%
  select(bug, resistance_drug, level)

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
dat = cons_res %>%
  rowwise %>%
  do(add_organisms(.)) %>%
  left_join(coherence, by=c('resistance_drug', 'bug')) %>%
  filter(level %in% c('state', 'hrr', 'zipcode'))

# do the correlations 