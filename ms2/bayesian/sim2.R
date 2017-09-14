#!/usr/bin/env Rscript

library(rstan)

inv_logit = function(x) 1/(1+exp(-x))

S = 30
Size = sample(1:5, S, replace=TRUE)
A = sum(Size)
abg_state = hashmap::hashmap(1:A, unlist(lapply(1:S, function(i) rep(i, Size[i]))))

cons = runif(S, min=0, max=0.5)

beta0 = 0.1
beta1 = 1
psi = 25
eta = beta0 + beta1 * cons
state_p = inv_logit(eta)
abg_p = lapply(1:S, function(s) rbeta(Size[s], psi*state_p[s], psi*(1 - state_p[s]))) %>%
  unlist

Iso = as.integer(rexp(A) * 1000)
Res = rbinom(A, Iso, abg_p)

dat = left_join(
  data_frame(state=1:S, cons, state_p),
  data_frame(abg=1:A, Iso, Res, state=abg_state[[abg]]), by='state') %>%
  mutate(empirical_p=Res/Iso)
# 
# n_bene = as.integer(rnorm(A, mean=1e5, sd=2e4))
# stopifnot(min(n_bene) > 0)
# 
# beta0 = 0.1
# beta1 = 1
# psi = 0.05
# 
# ineq = lapply(1:A, function(i) {
#   data_frame(n_claims=rpois(n_bene[i], cons[i])) %>%
#     count(n_claims) %>%
#     rename(n_bene=n) %>%
#     mutate(abg=i)
# }) %>% bind_rows %>%
#   group_by(abg) %>%
#   mutate(f=n_bene/sum(n_bene)) %>%
#   ungroup() %>%
#   mutate(true_p = beta0 + beta1 * n_claims)
# 
# abg = ineq %>%
#   group_by(abg) %>%
#   summarize(avg_p=weighted.mean(true_p, w=n_bene)) %>%
#   mutate(Iso=Iso) %>%
#   mutate(Res=rbinom(nrow(.), size=Iso, prob=avg_p))
# 
# NC = max(ineq$n_claims) + 1
# 
# Fmat = matrix(data=0, nrow=A, ncol=NC)
# 
# for (ineq_i in 1:nrow(ineq)) {
#   i = ineq$abg[ineq_i]
#   j = ineq$n_claims[ineq_i] + 1
#   f = ineq$f[ineq_i]
#   
#   Fmat[i, j] = f
# }
# 
# Res = abg$Res

fit = stan(file='model2.stan',
           data=c('S', 'A', 'Size', 'cons', 'Iso', 'Res'),
           iter=1000,
           chains=2)

pairs(fit, pars=c('lp__', 'psi', 'beta0', 'beta1'))