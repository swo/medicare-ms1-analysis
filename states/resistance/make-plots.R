#!/usr/bin/env Rscript

res = read_tsv('model-results-p005.tsv')
plot_dat = read_tsv('plotting-dat.tsv')

plot_f = function(this_bug, this_resistance_drug, this_explanatory, this_response) {
  this_res = res %>%
    filter(bug==this_bug, resistance_drug==this_resistance_drug, explanatory==this_explanatory, response==this_response)

  if (dim(this_res)[1] != 1) stop('wrong number of rows')

  z = this_res$correlation %>% head
  p.value = this_res$p.value %>% head

  plot_dat %>%
    filter(bug==this_bug, resistance_drug==this_resistance_drug) %>%
    ggplot(aes_string(this_explanatory, this_response)) +
    geom_point() +
    geom_smooth(method='lm', se=FALSE) +
    theme_minimal() +
    xlab(this_explanatory) +
    ylab(this_response) +
    ggtitle(sprintf('%s susceptibility to %s', this_bug, this_resistance_drug))
}

save_plot = function(bug, drug, expl, resp, file_nickname) {
  fn = sprintf('plots/%s.pdf', file_nickname)
  p = plot_f(bug, drug, expl, resp)
  ggsave(fn, plot=p, useDingbats=FALSE)
}

# dubious
save_plot('Coagulase-negative staphylococci', 'levofloxacin', 'claims_per_1k_ppl', 'min_sus', 'cnst') # why only min sus?
save_plot('Enterobacter cloacae', 'ceftriaxone', 'did', 'mean_sus', 'ecloacae') # why only DID?
save_plot('Escherichia coli', 'ampicillin', 'did', 'median_sus', 'ecoli-amp') # why only DID?
save_plot('Escherichia coli', 'trimethoprim-sulfamethoxazole', 'claims_per_1k_ppl', 'mean_sus', 'ecoli-tmpsmx') # why not claims per 1k ppl?
save_plot('Klebsiella pneumoniae', 'ceftriaxone', 'did', 'mean_sus', 'klebs') # why only DID?
save_plot('Methicillin-resistant S. aureus (MRSA)', 'levofloxacin', 'claims_per_1k_ppl', 'max_sus', 'mrsa-levo') # why only max sus?
save_plot('Methicillin-resistant S. aureus (MRSA)', 'trimethoprim-sulfamethoxazole', 'did', 'mean_sus', 'mrsa-tmpsmx') # why only DID?
save_plot('Proteus mirabilis', 'levofloxacin', 'claims_per_1k_ppl', 'min_sus', 'proteus') # why only min sus?

# middling
save_plot('Streptococcus pneumoniae (non-meningitis)', 'trimethoprim-sulfamethoxazole', 'claims_per_1k_ppl', 'mean_sus', 'pneumo-tmpsmx') # why only some?

# likely
save_plot('Escherichia coli', 'ciprofloxacin', 'claims_per_1k_ppl', 'mean_sus', 'ecoli-cipro')
save_plot('Escherichia coli', 'levofloxacin', 'claims_per_1k_ppl', 'mean_sus', 'ecoli-levo')
save_plot('Methicillin-susceptible S. aureus (MSSA)', 'ciprofloxacin', 'claims_per_1k_ppl', 'mean_sus', 'mssa-cipro')
save_plot('Methicillin-susceptible S. aureus (MSSA)', 'erythromycin', 'claims_per_1k_ppl', 'mean_sus', 'mssa-aziery')
save_plot('Streptococcus pneumoniae (non-meningitis)', 'erythromycin', 'claims_per_1k_ppl', 'mean_sus', 'pneumo-aziery')
