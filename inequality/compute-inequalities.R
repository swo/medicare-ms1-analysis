# summarize unnevenness

library(ineq)
library(VGAM)

olesen_index = function(counts) {
  x = data.frame(counts=counts) %>%
    arrange(counts) %>%
    mutate(x=1:dim(.)[1], x=x/max(x), y=cumsum(counts)/sum(counts)) %>%
    filter(x + y - 1 < 0) %$%
    x
  
  if (length(x) > 0) {
    tail(x, 1)
  } else {
    1.0
  }
}

printf = function(fmt, ...) sprintf(fmt, ...) %>% cat
gini = function(x) ineq(x, type='Gini')
user_gini = function(x) gini(x[x > 0])

power_law_fit = function(counts) {
  x = counts[counts > 0]
  const = 1 / length(x) * sum(log(x))
  g = function(a) zeta(a, 1) / zeta(a) + const
  uniroot(g, c(1 + 1e-6, 100))$root
}

add_zeros = function(n_total, counts) {
  nonzero_counts = counts[counts > 0]
  n_zero_benes = n_total_benes - length(nonzero_counts)
  c(nonzero_counts, rep(0, times=n_zero_benes)) %>% .[order(.)]
}

n_total_benes = 10341351

plot_cdf = function(counts) {
  data.frame(counts=counts) %>%
    arrange(counts) %>%
    mutate(x=1:dim(.)[1], x=x/max(x), y=cumsum(counts)/sum(counts)) %>%
    sample_frac(1e-3) %>%
    ggplot(aes(x=x, y=y)) + geom_line() + coord_fixed()
}

ineq_table = function(wap) {
  abxs = wap %>% select(-bene) %>% names
  
  gini = function(x) ineq(x, type='Gini')
  
  n_abxs = length(abxs)
  users = rep(NA, length=n_abxs)
  days = rep(NA, length=n_abxs)
  ginis = rep(NA, length=n_abxs)
  
  for (i in 1:n_abxs) {
    abx = abxs[i]
    counts = wap[[abx]] %>% .[. > 0]
    
    users[i] = length(counts)
    days[i] = sum(counts)
    
    if (users[i] > 1) {
      ginis[i]= gini(counts)
      alphas[i] = power_law_fit(counts)
    }
  }
  
  res = data.frame(abx=abxs, users=users, days=days, gini=ginis, alpha=alphas)
  res
}