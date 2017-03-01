# Source the state-report.Rmd before doing any of this

# There are a few ways to do coherence:
# - Kruskal-Wallis
# - ANOVA on a linear model
# - Permutation tests using F statistics or SSE (my previous approach)
# - Dumb things like number of antibiograms or states

# Kruskal Wallis

safe_kruskal = function(x) {
  try(return(kruskal.test(percent_nonsusceptible ~ factor(state), data=x) %>% tidy), silent=TRUE)
  data_frame(statistic=NA, p.value=1.0, parameter=NA, method=NA)
}

state_dat %>%
  group_by(bug, drug_group) %>%
  do(safe_kruskal(.)) %>%
  mutate(significant=p.value < 0.05) %>%
  ggplot(aes(x=p.value, fill=significant)) +
  geom_histogram() +
  scale_x_log10()

# ANOVA

lm_anova = function(x) {
  try(return(
    lm(percent_nonsusceptible ~ factor(state), data=x) %>%
      anova %>% tidy %>%
      filter(term=='factor(state)') %>%
      select(-term)
  ), silent=TRUE)
  data_frame(df=NA, sumsq=NA, meansq=NA, statistic=NA, p.value=1.0)
}

state_dat %>%
  group_by(bug, drug_group) %>%
  do(lm_anova(.)) %>%
  mutate(significant=p.value < 0.05) %>%
  ggplot(aes(x=p.value, fill=significant)) +
  geom_histogram() +
  scale_x_log10()

# Permutation tests

F_statistic = function(values, groups) {
  grand_mean = mean(values)
  group_means = data_frame(value=values, group=groups) %>%
    group_by(group) %>%
    summarize(group_mean=mean(value)) %$%
    group_mean
  
  # sum of squares of groups
  ssg = sum((group_means - grand_mean) ** 2)
  sse = sse_statistic(values, groups)
  
  n_groups = length(unique(groups))
  n_values = length(values)
  (ssg / (n_groups - 1)) / (sse / (n_values - 1))
}

shuffled = function(x) sample(xtfrm(x))
perm_test = function(values, groups, statistic_f, n_permutations=1e3) {
  n_smaller = 0
  true_statistic = statistic_f(values, groups)
  for (i in 1:n_permutations) {
    this_statistic = statistic_f(values, shuffled(groups))
    if (this_statistic < true_statistic) n_smaller = n_smaller + 1
  }
  (n_smaller + 1) / (n_permutations + 1)
}

sse_statistic = function(values, groups) {
  data_frame(value=values, group=groups) %>%
    group_by(group) %>%
    mutate(group_mean=mean(value)) %>%
    ungroup %>%
    mutate(within_group_residual=value-group_mean) %$%
    sum(within_group_residual ** 2)
}

z = state_dat %>%
  group_by(bug, drug_group) %>%
  summarize(p.value=perm_test(percent_nonsusceptible, state, sse_statistic, n_permutations=1e2)) %>%
  ungroup()

z %>%
  mutate(significant=p.value < 0.05) %>%
  ggplot(aes(x=p.value, fill=significant)) +
  geom_histogram() +
  scale_x_log10()

# Comparing methods

x = state_dat %>% group_by(bug, drug_group) %>% count(state) %>% summarize(n_abg=sum(n), n_states=n())
y = state_dat %>% group_by(bug, drug_group) %>% do(safe_kruskal(.))
z = state_dat %>% group_by(bug, drug_group) %>% do(lm_anova(.))
a = left_join(x, y, by=c('bug', 'drug_group')) %>%
  left_join(z, by=c('bug', 'drug_group'), suffix=c('.kw', '.an'))