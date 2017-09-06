data {
  int<lower=0> S; // number of states
  int<lower=0> A; // number of antibiograms (total)
  int<lower=0> Size[S]; // number of antibiograms in each state
  int<lower=0> Iso[A]; // number of isolates in each antibiogram
  int<lower=0> Res[A]; // number of resistant isolates in each antibiogram
  vector<lower=0>[S] Cons; // consumption
}

parameters {
  vector<lower=0, upper=1>[A] p; // probability of resistance in each antibiogram
  vector<lower=0, upper=1>[S] mu; // mean resistance in each state
  vector<lower=0>[S] phi; // narrowness of resistance in each state
  real beta0; // "intercept" in logistic function
  real beta1; // "slope" in logistic function
  real<lower=0> psi; // narrowness of the consumption-resistance relationship
}

transformed parameters {
  vector[S] eta; // linear predictor
  vector<lower=0, upper=1>[S] muhat; // estimator
  eta = beta0 + beta1 * Cons;
  muhat = inv_logit(eta);
}

model {
  // initialize counter to point to position in antibiogram data vectors while
  // looping over states
  int pos;
  pos = 1;

  // For the antibiogram data and parameters (Iso, Res, p), we go state by state
  // because the data is ragged (different number of abg in each state)
  // cf. p. 231 in Stan reference for using this construct
  for (state in 1:S) {
    segment(Res, pos, Size[state]) ~ binomial(segment(Iso, pos, Size[state]), p[state]);
    // In the beta distributions, alpha=mu*phi, beta=(1-mu)*phi
    segment(p, pos, Size[state]) ~ beta(phi[state] * mu[state], phi[state] * (1-mu[state]));
    pos = pos + Size[state];
  }

  phi ~ cauchy(0, 25); // mean, sigma
  beta0 ~ normal(0, 10);
  beta1 ~ normal(0, 10);
  mu ~ beta(psi * muhat, psi * (1-muhat));
  psi ~ cauchy(0, 50);
}
