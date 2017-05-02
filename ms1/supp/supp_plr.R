plr = function(test, disease) {
  # positive likelihood ratio:
  # P(T+ | D+) / P(T+ | D-)
  # = N_true_pos * N_false / (N_false_pos * N_true)
  
  nsum = function(x) as.numeric(sum(x)) # to prevent overflow
  n_true_pos = nsum(test & disease)
  n_false_pos = nsum(test & !disease)
  n_true = nsum(disease)
  n_false = nsum(!disease)
  
  #n_true_pos * n_true / (n_false_pos * n_false)
  n_true_pos * n_false / (n_false_pos * n_true)
}

plr_matrix = function(df) {
  crossing(test_condition=names(df), true_condition=names(df)) %>%
    filter(test_condition != true_condition) %>%
    rowwise %>%
    mutate(PLR=plr(df[[test_condition]], df[[true_condition]]))
}

condition_names = read_tsv('../../chronic-conditions/names.txt',
                           col_names=c('condition', 'condition_name'))

# load data, keeping only the CC data for the qualifying beneficiaries
bene_ids = read_tsv('../../bene_2011.tsv') %>%
  filter(between(age, 67, 96), hmo_months==0) %$%
  bene_id

cc = read_feather('../../cc_2011.feather') %>%
  filter(bene_id %in% bene_ids) %>%
  select(-ALZH)

m = cc %>% select(-bene_id) %>% plr_matrix

conditions = c("AMI", "ISCHMCHT", "CHF", "HYPERT", "CHRNKIDN", "STRKETIA", "HYPERL", 
               "COPD", "ASTHMA", "CNCRLUNG",
               "HIPFRAC", "OSTEOPRS", "ANEMIA", "DEPRESSN", 
               "ALZHDMTA", "DIABETES", "ATRIALFB", "RA_OA", 
               "HYPOTH", "CNCRBRST",  "CNCRENDM",  
               "CNCRCLRC", "CNCRPRST", "HYPERP",
               "CATARACT", "GLAUCOMA")
male_conditions = c('CNCRPRST', 'HYPERP')
female_conditions = c('CNCRENDM', 'CNCRBRST', 'OSTEOPRS', 'HYPOTH')

plot_f = function(cc, excluded_conditions) {
  levels = setdiff(conditions, excluded_conditions)
  condition_factor = function(x) factor(x, levels=levels, ordered=TRUE)
  
  cc %>%
    filter(test_condition %in% levels, true_condition %in% levels) %>%
    mutate_at(vars(test_condition, true_condition), condition_factor) %>%
    ggplot(aes(x=true_condition, y=test_condition, fill=log(PLR))) +
    geom_tile() +
    scale_fill_gradient2() +
    scale_y_discrete(limits=rev(levels)) +
    xlab('predicted condition') + ylab('predictor condition') +
    theme(axis.text.x=element_text(angle=90, hjust=1),
          panel.background=element_rect(fill='black'),
          panel.grid=element_blank())
}

male_plot = plot_f(m, female_conditions)
ggsave('plr_male.pdf', plot=male_plot, useDingbats=FALSE)
female_plot = plot_f(m, male_conditions)
ggsave('plr_female.pdf', plot=female_plot, useDingbats=FALSE)