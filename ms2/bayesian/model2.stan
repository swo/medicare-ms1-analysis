data {
  int<lower=0> S;
  int<lower=0> A;
  int<lower=0> NC;
  int<lower=0> Size[S];
  vector<lower=0>[S] cons;
  int<lower=0> Iso[A];
  int<lower=0> Res[A];
  matrix<lower=0>[S, NC] Fmat;
}

parameters {
  real beta0; // regression intercept
  real beta1; // regression coefficient
  real<lower=0> psi;
  vector<lower=0, upper=1>[A] p;
}

transformed parameters {
  vector[NC] eta;
  vector[S] rho;
  
  for (n in 1:NC) {
    eta[n] = inv_logit(beta0 + beta1 * n);
  }
  
  rho = Fmat * eta;
}

model {
  int pos;
  pos = 1;
  
  for (s in 1:S) {
    segment(p, pos, Size[s]) ~ beta(psi*rho[s], psi*(1-rho[s]));
    pos = pos + Size[s];
  }
  
  Res ~ binomial(Iso, p);
  
  // priors
  beta0 ~ normal(0, 5);
  beta1 ~ normal(0, 5);
  psi ~ normal(0, 5);
}
