data {
  int<lower=0> S; // number of states
  int<lower=0> A; // number of antibiograms (total)
  int<lower=0> Size[S]; // number of antibiograms in each state
  int<lower=0> Iso[A]; // number of isolates in each antibiogram
  int<lower=0> Res[A]; // number of resistant isolates in each antibiogram
  vector<lower=0>[S] Fnz; // fraction nonzero
  vector<lower=0>[S] Mup; // mean consumption among consumers
}

parameters {
  vector<lower=0, upper=1>[A] p; // probability of resistance in each abg
  real beta0; // regression intercept
  real beta_fnz; // regression coefficient
  real beta_mup;
  real<lower=0> psi; // narrowness of beta predictive distribution
}

transformed parameters {
  vector[S] theta; // linear predictor
  theta = beta0 + beta_fnz * Fnz + beta_mup * Mup;
}

model {
  int pos;
  pos = 1;

  // For the antibiogram data and parameters (Iso, Res, p), we go state by state
  // because the data is ragged (different number of abg in each state)
  // cf. p. 231 in Stan reference for using this construct
  for (s in 1:S) {
    // The p's for abgs in a state are down from a beta distribution whose mean
    // is predicted by the consumption data (and whose variance is variable)
    segment(p, pos, Size[s]) ~ beta(theta[s] * psi, (1-theta[s]) * psi);
    pos = pos + Size[s];
  }

  Res ~ binomial(Iso, p); // # resistants is determined by the p's
  beta0 ~ normal(0, 5); // weak priors on the regression coefficients
  beta_fnz ~ normal(0, 5);
  beta_mup ~ normal(0, 5);
  psi ~ cauchy(0, 5); // weak prior on the narrowness ("strength") of predictive beta
}
