data {
  int<lower=0> S; // number of states
  int<lower=0> A; // number of antibiograms (total)
  int<lower=0> Size[S]; // number of antibiograms in each state
  int<lower=0> Iso[A]; // number of isolates in each antibiogram
  int<lower=0> Res[A]; // number of resistant isolates in each antibiogram
  vector<lower=0>[S] Cons; // consumption
}

parameters {
  vector<lower=0>[A] p; // probability of resistance in each antibiogram
  vector<lower=0>[S] mu; // mean resistance in each state
  vector<lower=0>[S] phi; // narrowness of resistance in each state
  real<lower=0> g0; // consumption leading to 50% resistance
  real<lower=0> g1; // 1/g1 is the fraction above g0 leading to 75% resistance
  real<lower=0> psi; // narrowness of the consumption-resistance relationship
}

transformed parameters {
  vector[S] eta; // linear predictor
  vector[S] muhat; // estimator
  eta = g1 * ((Cons / g0) - 1);
  muhat = 1.0 ./ (1 + exp(-eta));
}

model {
  // initialize counter to point to position in antibiogram data vectors while
  // looping over states
  int pos;
  pos = 1;

  // cf. p. 231 in Stan reference for using this construct
  for (state in 1:S) {
    segment(Res, pos, Size[state]) ~ binomial(segment(Iso, pos, Size[state]), p[state]);
    pos = pos + Size[state];
  }

  p ~ beta(phi .* mu, phi .* (1-mu)); // alpha=mu*phi, beta=(1-mu)*phi
  phi ~ cauchy(0, 25); // mean, sigma
  g0 ~ cauchy(0, 50);
  g1 ~ cauchy(0, 10);
  mu ~ beta(psi * muhat, psi * (1-muhat));
  psi ~ cauchy(0, 50);
}
