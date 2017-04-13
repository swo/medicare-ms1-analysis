# > names(x)
# [1] "state"    "cpb"      "v"        "fnznc"    "nzgini"   "n_bene"   "n_claims" "fz"      

x = usage %>%
  group_by(state) %>%
  summarize(cpb=mean(n_claims),
            v=var(n_claims),
            fnznc=fnz(n_claims),
            fz=1-fnznc,
            nzgini=Gini(nz(n_claims)),
            n_bene=n(),
            n_claims=sum(n_claims))

# fit the p value
m = nls(fz ~ (1-p) ** ((1-p)/p*cpb), data=x, start=list(p=0.5))
p = coef(m)[['p']]

# or, do this with the variance
# the coefficient of the variance v is (1-p)
# this give 0.251 -> 0.749, which is higher than p=0.7 from above
# lm(cpb ~ v + 0, data=x)

# get the r values
r = (1-p)/p * x$cpb

gini_range = function(r, n_individuals, n_trials=1e2) {
  values = replicate(n_trials, Gini(nz(rnbinom(n_individuals, r, p))))
  data_frame(mu=mean(values), cil=quantile(values, 0.05), ciu=quantile(values, 0.95))
}

y = data_frame(r=r, n_bene=x$n_bene) %>%
  group_by(r) %>%
  do(gini_range(.$r, .$n_bene)) %>%
  mutate(nbmu=p*r/(1-p), nbf0=(1-p)**r)