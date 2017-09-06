# Data

- \# states $S$
- \# antibiograms in each state $A_s$
- \# isolates in each abg $I_{sa}$
- \# nonsusceptible isolates in each abg $R_{sa}$
- consumption in each state $C_s$
- (optional) variance in consumption in each state

# Parameters

- probability of resistance $p_{sa}$ for each abg
- mean $\mu_s$ for each state
- measure $\phi_s$ ($\equiv \alpha_s + \beta_s$) of the variance in the resistances for each state
- $\gamma_0$, $\gamma_1$ overall
- $\psi$, a measure of the variance in the consumption-resistance relationship

# Transformed parameters

I parameterize the logit a little differently. Typically, we write
$\mathrm{logit}(\beta_0 + \beta_1 x)$ such that
$\beta_1 x = -\beta_0 \implies y = \tfrac{1}{2}$
and $\beta_1 x = -\beta_0 + 1 \implies y \approx \tfrac{3}{4}$. Instead, I write
$\mathrm{logit}[\gamma_1 (x/\gamma_0 - 1)]$ such that
$x = \gamma_0 \implies y=\tfrac{1}{2}$ and
$x = (1 + 1/\gamma_1) \gamma_0 \implies y\approx\tfrac{3}{4}$.
Thus, $\gamma_0$ is the 50% midpoint and $1/\gamma_1$ is the fraction
above $\gamma_0$ that you have to go to get to 75%.
(It's actually $e/(e+1) = 0.73$.)

In the priors, we will prefer small $\gamma_0$ (i.e., we think that the
consumption midpoint is not arbitrarily high) and small $\gamma_1$ (i.e., we
think that the slope at the midpoint is not arbitrarily high).

- linear predictor $\eta_s = \gamma_1 (C_s/\gamma_0 - 1)$
- estimator $\hat{\mu}_s = \mathrm{logit}(\eta_s)$

# Model

Under the parameterization of the beta distribution that we're using, namely
$\mathrm{Beta}(\mu, \phi)$, the variance is $\mu (1-\mu) / (\phi + 1)$. The variance
is maximized for $\mu = \tfrac{1}{2}$, so that the variance is at most
$1/[4(\phi + 1)]$. For the scale parameters $\phi_s$ and $\psi$, we set the
prior using a Half-Cauchy distribution. The hyperparameter in the Half-Cauchy is
the position at which the cdf is $\tfrac{1}{2}$, i.e., it is as likely that the
value will be below as above this position.

I'm not sure what to use for the prior on the central position $\gamma_0$. Maybe
Half-Cauchy isn't the right thing.

- $R_{sa} \sim \mathrm{Binom}(p_{sa}, I_{sa})$ (i.e., resistants are drawn as yes-no trials)
- $p_{sa} \sim \mathrm{Beta}(\mu_s, \phi_s)$ (i.e., yes-no probability for different hospitals is drawn from a beta distribution)
- $\phi_s \sim \mathrm{HalfCauchy}(\tau_\phi)$ (i.e., the priors on the narrowness of the state-level resistance distributions)
- $\gamma_0 \sim \mathrm{HalfCauchy}(\tau_{\gamma_0})$ (prior on the central position of the consumption-resistance relationship)
- $\gamma_1 \sim \mathrm{HalfCauchy}(\tau_{\gamma_1})$ (prior on the scale)
- $\mu_s \sim \mathrm{Beta}(\hat{\mu}_s, \psi)$ (i.e., we expect to draw the state-level resistances from a beta distribution whose mean is the thing predicted by the logistic model and whose variance comes from an overall spread parameter $\psi$)
- $\psi \sim \mathrm{HalfCauchy}(\tau_\psi)$ (prior on the variance in the consumption-resistance relationship)

# Values for the hyperparameters

- $\tau_\phi$: I expect a standard deviation on the order of 10% within a state. That means a variance of $10^{-2}$, which would be $\phi = 25$. So let's put this at 25.
- $\tau_{\gamma_0}$: This could depend strongly on the drug. Instead let's say that if consumption were increased 2-fold from the highest-consuming state, would we expect 50% resistance? Maybe. 10-fold? Seems very likely. So let's put this at 10-fold higher than the highest-consumed drug, which will be an overestimate for other drugs. Mean consumption for quinolones tops out at 0.5 claims per beneficiary, so let's put this at 5.
- $\tau_{\gamma_1}$: If you consumption needs to be $x$ to get 50% resistance, what does consumption have to be to get 75% resistance? Maybe $2x$, but it seems unlikely that it's $100x$. So let's put this at $10$.
- $\tau_\psi$: How much variance in the consumption-resistance relationship do we expect? If you had multiple states with exactly the consumption that will give rise to a mean 50% resistance, what variation do we expect? I think this could be a little higher than for $\phi$, so let's also put this at 50.
