# True distribution in each state

Each state has some true distribution of susceptibilities (or resistances). Let
each hospital $i$ in the state have probability (of, say, susceptibility) $p_i$.
Then, given the "successes" $s_i$ and of $n_i$ total "trials" (i.e., susceptible
isolates of all isolates) for each hospital, we aim to infer the distribution in
the states, treating the actual probabilities $p_i$ as nuisance parameters.

The likelihood of the data $X \equiv \{s_i, n_i\}_i$ given the parameters
$\alpha, \beta$ of the true beta distribution within the state is:

$$
\begin{aligned}
P(X | \alpha, \beta) &\propto
  \prod_i \int_0^1 \mathrm{Bin}(s_i; n_i, p_i) \mathrm{Beta}(p_i; \alpha, \beta) \,\mathrm{d}p_i \\
&= \prod_i \int_0^1 \binom{n_i}{s_i} p_i^{s_i} (1-p_i)^{n_i-s_i}
  \frac{1}{B(\alpha, \beta)} p_i^{\alpha - 1} (1-p_i)^{\beta-1} \,\mathrm{d}p_i \\
&= \prod_i \binom{n_i}{s_i} \frac{1}{B(\alpha, \beta)}
  \int_0^1 p_i^{s_i+\alpha-1} (1-p_i)^{n_i-s_i+\beta-1} \,\mathrm{d}p_i \\
&= \prod_i \binom{n_i}{s_i} \frac{1}{B(\alpha, \beta)} B(s_i + \alpha, n_i - s_i + \beta)
\end{aligned}
$$

where $B$ is the beta function. The posterior probability of the parameters
$\alpha$ and $\beta$ is

$$
P(\alpha, \beta | X) \propto P(\alpha, \beta) P(X | \alpha, \beta)
$$

If we assume a flat, improper prior $P(\alpha, \beta) = 1$, then the maximum
likelihood estimates for $\alpha$ and $\beta$ are the ones that maximize

$$
f(\alpha, \beta) = \prod_i \frac{B(s_i + \alpha, n_i - s_i + \beta)}{B(\alpha, \beta)}
$$

or, equivalently, the logarithm

$$
\log f(\alpha, \beta) = \sum_i \log B(s_i + \alpha, n_i - s_i + \beta) - N \log B(\alpha, \beta)
$$

where $N$ is the number of hospitals in the state.

# Inference of the correlation coefficient

## Model setup

The real question is: given the susceptible and resistance isolates from each
hospital in each state, and given the consumption in each state, how does
consumption affect resistance?

So now we have a three-level model:

1. Consumption determines the state's true hyper-mean resistance
2. The state's distribution on resistances determines the hospital's resistances
3. The hospital's resistances determine the outcomes of the isolate Bernoulli trials

The number of isolates and the consumption per state are fixed values from the data.
The state's true hyper-mean variance and the individual hospital resistances are
treated as nuisance parameters.

This approach involves a re-parameterization: rather than considering $\alpha$ and
$\beta$ for the state-level distribution, we are considering the distribution's mean
$\mu$ and variance $\sigma$. By reversing the well-known equations

$$
\begin{gathered}
\mu = \frac{\alpha}{\alpha + \beta} \\
\sigma = \frac{\alpha \beta}{(\alpha + \beta)^2 (\alpha + \beta + 1)}
\end{gathered}
$$

we obtain

$$
\begin{aligned}
\alpha &= \frac{1}{\sigma}(\mu^2 - \mu^3 - \mu\sigma) \\
\beta &= \frac{1}{\sigma}(\mu - 1)(\mu^2 - \mu + \sigma)
\end{aligned}
$$

Similarly, the parameters $k$ and $\theta$ in a gamma distribution can
be related to the distribution's mean $\nu$ and variance $\tau$:

$$
\begin{gathered}
\nu = k \theta \\
\tau = k \theta^2
\end{gathered}
$$

from which it follows that

$$
\begin{gathered}
k = \tau / \nu^2 \\
\theta = \tau / \nu
\end{gathered}
$$

Let $C$ be the consumption in a state, $m$ be the "slope" of the
consumption-resistance line, $p$ be the probability of resistance in that
state, $\mu$ and $\sigma$ be the mean and variance of the true resistance
distribution, $\nu$ and $\tau$ be the mean and variance of the true
distribution of state-level means. Then the probability of the data in a
hospital (really just the $s$) is

$$
P(X | m, p, \mu, \sigma, \tau) =
  \mathrm{Bin}(s; n, p) \mathrm{Beta}(p; \mu, \sigma) \mathrm{Gamma}(\mu; \nu, \tau)
$$

where $\nu = mC$.


## Inference

We treat the consumption $C$ and number of isolates $n_i$ in each hospital as
fixed. We intend to only do inference on $m$, and we treat the other variables
$p_i$, $\mu$, $\sigma$, and $\tau$ as nuisances. Thus, for a particular state,
assuming a flat prior,

$$
P(m | X) = \prod_i \iiiint
  \mathrm{Bin}(s_i; n_i, p_i) \mathrm{Beta}(p_i; \mu, \sigma) \mathrm{Gamma}(\mu; Cm, \tau)
  \,dp_i \,d\mu \,d\sigma \,d\tau
$$

For the full scenario, we need an extra set of indices over states $j$:

$$
\begin{split}
P(m | X) = \prod_{i,j} \iiiint
  \mathrm{Bin}(s_{i,j}; n_{i,j}, p_{i,j}) \mathrm{Beta}(p_{i,j}; \mu_j, \sigma_j) \times \\
  \mathrm{Gamma}(\mu_j; C_j m, \tau_j) \,dp_{i,j} \,d\mu_j \,d\sigma_j \,d\tau_j
\end{split}
$$

This presents a difficulty: the integrals will be very small. It's probably
possible to do.

# Another parameterization

We assume that consumption correlates with resistance according to $\mu_i = m C_i$,
where $\mu_i$ is the average resistance in state $i$, $C_i$ is the
resistance in that state, and $m$ is the slope of the linear relationship.

We assume that each state has some true distribution of hospital-level
resistances $p_{ij}$, which are nuisance parameters. Those resistance
probabilities are drawn from a beta distribution with mean $m C_i$ and a
variance $V_i \in (0, \tfrac{1}{4}]$ that we also consider a nuisance
parameter.

Each hospital $j$ in state $i$ has $n_{ij}$ isolates for a bug/drug combination.
The probability of finding $s_{ij}$ "successes" among those $n_{ij}$ isolates is
determined by the binomial distribution.

Thus, the likelihood of the data, the $s_{ij}$, given the undetermined parameter
of interest $m$ is

$$
P(\{s_{ij}\} | m) = \prod_i \int_{V_i} \prod_j \int_{p_{ij}}
  \mathrm{Beta}(p_{ij}; \alpha_i, \beta_i) \mathrm{Bin}(s_{ij}; n_{ij}, p_{ij})
  \,\mathrm{d}p_{ij} \,\mathrm{d}V_i
$$

where the parameters in the beta function are

$$
\begin{gathered}
\alpha_i = \frac{1}{V_i}(\mu_i^2 - \mu_i^3 - \mu_i V_i) \\
\beta_i = \frac{1}{V_i}(\mu_i - 1)(\mu_i^2 - \mu_i + V_i)
\end{gathered}
$$

Note that

$$
\begin{gathered}
\mathrm{Beta}(p; \alpha, \beta) = \frac{1}{\mathrm{B}(\alpha, \beta)}
  p^{\alpha-1} (1-p)^{\beta-1} \\
\mathrm{Bin}(s; n, p) = \binom{n}{s} p^s (1-p)^{n-s}
\end{gathered}
$$

where $\mathrm{B}(\alpha, \beta)$ is the beta function, from which
we see that we can combine the powers of $p$ and $1-p$:

$$
P(\{s_{ij}\} | m) \propto \prod_i \int_{V_i}
  [\mathrm{B}(\alpha_i, \beta_i)]^{-N_i}
  \prod_j \mathrm{B}(\alpha_i + s_{ij}, \beta_i + n_{ij} - s_{ij})
  \,\mathrm{d}V_i
$$

where $N_i$ is the number of hospitals in state $i$.

We can remove the first product by taking the logarithm:

$$
\log P(\{s_{ij}\} | m) \propto \sum_i \log\left\{ \int_{V_i}
  [\mathrm{B}(\alpha_i, \beta_i)]^{-N_i}
  \prod_j \mathrm{B}(\alpha_i + s_{ij}, \beta_i + n_{ij} - s_{ij})
  \,\mathrm{d}V_i \right\}
$$
