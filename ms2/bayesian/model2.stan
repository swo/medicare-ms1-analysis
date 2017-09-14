data {
  int<lower=0> A; // number of antibiograms (total)
  vector<lower=0>[A] cons;
  int<lower=0> Iso[A];
  int<lower=0> Res[A];
}

parameters {
  real beta0; // regression intercept
  real beta1; // regression coefficient
  real<lower=0> psi;
  vector<lower=0, upper=1>[A] p;
}

transformed parameters {
  vector[A] rho;
  rho = beta0 + beta1 * cons;
}

model {
  Res ~ binomial(Iso, rho);
  
  // priors
  beta0 ~ normal(0, 5);
  beta1 ~ normal(0, 5);
  psi ~ cauchy(0, 5);
}
