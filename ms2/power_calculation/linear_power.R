#!/usr/bin/env Rscript

f = function(n, b1, x1.mu, x1.sd, b2, x2.mu, x2.sd, e.sd) {
  n = 40
  x1.mu = 0.082
  x2.mu = 0.309
  x1.sd = 0.016
  x2.sd = 0.027
  e.sd = 4.96

  b1 = 115.6
  b2 = -99.77
  b0 = 0.0

  x1 = rnorm(n, sd=x1.sd)
  x2 = rnorm(n, sd=x2.sd)
  e = rnorm(n, sd=e.sd)

  y = b0 + b1 * x1 + b2 * x2 + e
  m = lm(y ~ x1 + x2)
  return(m %>% tidy %>% filter(term=='x2') %$% p.value)
}

power = function(n_trials, p.value.f, alpha=0.05) {
  p.values = replicate(n_trials, p.value.f())
  sum(p.values < alpha) / n_trials
}

ecoli_bactrim_f = function() f(44, 115.6, 0.082, 0.016, -99.77, 0.309, 0.027, 4.96)
ecoli_betalactam_f = function() f(44, 269.74, 0.158, 0.015, -187.42, 0.270, 0.015, 8.6)

show(power(1e3, ecoli_bactrim_f))
show(power(1e3, ecoli_betalactam_f))
