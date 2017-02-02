#!/usr/bin/env Rscript

res = read_tsv('model-results-p005.tsv')
plot_dat = read_tsv('plotting-dat.tsv')

plot_f = function(this_consumption_drug, this_resistance_drug, this_explanatory, this_response) {
  this_res = res %>%
    filter(consumption_drug==this_consumption_drug, resistance_drug==this_resistance_drug, explanatory==this_explanatory, response==this_response)

  if (dim(this_res)[1] != 1) stop('wrong number of rows')

  z = this_res$correlation %>% head
  p.value = this_res$p.value %>% head

  plot_dat %>%
    filter(consumption_drug==this_consumption_drug, resistance_drug==this_resistance_drug) %>%
    ggplot(aes_string(this_explanatory, this_response)) +
    geom_point() +
    geom_smooth(method='lm', se=FALSE) +
    theme_minimal() +
    xlab(this_explanatory) +
    ylab(this_response) +
    ggtitle(sprintf('Pneumo. sus. to %s, cons. of %s', this_resistance_drug, this_consumption_drug))
}

save_plot = function(cons, res, expl, resp, file_nickname) {
  fn = sprintf('plots/%s.pdf', file_nickname)
  p = plot_f(cons, res, expl, resp)
  ggsave(fn, plot=p, useDingbats=FALSE)
}

save_plot('penicillin', 'penicillin', 'claims_per_1k_ppl', 'mean_sus', 'pen')
