import::from(gplots, heatmap.2)

load('correlations.RData')

# Heatmap. Sqrt of MCC used. Female conditions are excluded.
female_conditions = c('CNCRENDM', 'CNCRBRST', 'OSTEOPRS')
pdf('mcc.pdf')
(!(rownames(correlations) %in% female_conditions)) %>%
  { correlations[., .] } %>%
  apply(1, function(x) pmax(x, 0)) %>%
  sqrt %>%
  gplots::heatmap.2(trace='none', symm=TRUE)
dev.off()
