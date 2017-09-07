# Data

- no. states $S$
- no. antibiograms in each state $A_s$
- no. isolates in each abg $I_{sa}$
- no. nonsusceptible isolates in each abg $R_{sa}$
- consumption in each state $C_s$ (or, $f_>$ and $\mu_>$)
- (*maybe later*) variance in consumption in each state

# Parameters

- probability of resistance $p_{sa}$ for each abg
- $\beta_0$, $\beta_1$: regression coefficients
- $\psi$, a measure of the variance in the consumption-resistance relationship

# Transformed parameters

- linear predictor $\eta_s = \beta_0 + \beta_1 C_s$
- estimator $\theta$: Now I'm just doing a linear regression, so $\theta_s = \eta_s$

# Model

In essence, this is a beta regression with an identity link function: the
estimator is just the linear predictor, but the errors on the estimates are
beta-distributed.

- The antibiogram data inform the resistance probabilities $p$ via: $R_{sa} \sim \mathrm{Binom}(p_{sa}, I_{sa})$ (i.e., resistants are drawn as yes-no trials)
- The resistance probabilities and consumption data inform the regression coefficients via: $p_{sa} \sim \mathrm{Beta}\left[\theta_s \psi, (1-\theta_s) \psi \right]$
- There are weak priors on the regression parameters:
    - $\beta_0 \sim \mathcal{N}(0, 5)$
    - $\beta_1 \sim \mathcal{N}(0, 5)$
    - $\psi \sim \mathrm{Half-}\mathcal{N}(0, 5)$
