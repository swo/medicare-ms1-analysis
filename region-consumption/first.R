abx = read_tsv("~/grad/proj/medicare/analysis/bene-abx-2011.tsv")
census = read_tsv("~/grad/proj/medicare/db/census-regions/census-regions.tsv")

# merge in the census region data
abx %<>%
  rename(state_abbreviation=state) %>% 
  left_join(census %>% select(-state))

counts_table = function(x) {
  counts = tabulate(x)
  n = 1:length(counts)
  data_frame(n=n, counts=counts)
}

plot_days = function(days) {
  days %>%
    counts_table %>% 
    ggplot(aes(x=n, y=counts)) + geom_point() + scale_y_log10() + xlab("# days") + ylab("# beneficiaries")
}

plot_drug_region = function(abx, drug, region) {
  abx %>% 
    rename_(.dots=setName(drug, 'drug')) %>%
    filter(region==region) %$%
    plot_days(drug)
}

plot_drug_region2 = function(abx, drug, region1, region2) {
  a = abx %>% 
    rename_(.dots=setNames(drug, 'drug')) %>% 
    filter(region==region1) %$%
    counts_table(drug)
  
  b = abx %>% 
    rename_(.dots=setNames(drug, 'drug')) %>% 
    filter(region==region2) %$%
    counts_table(drug)
  
  color1 = 'black'
  color2 = 'blue'
  title = sprintf('%s consumption in %s (%s) and %s (%s)', drug, region1, color1, region2, color2)
  ggplot(a, aes(x=n, y=counts)) +
    geom_point() + geom_smooth(col='black', se=FALSE) +
    geom_point(dat=b, col='blue') + geom_smooth(dat=b, col='blue', se=FALSE) +
    scale_y_log10() +
    xlab('# days') + ylab('# beneficiaries') + ggtitle(title) +
    theme_minimal()
}
