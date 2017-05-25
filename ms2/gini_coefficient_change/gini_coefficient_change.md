---
title: 'Changes in Gini coefficient and fraction nonzero'
author: 'Scott Olesen'
header-includes:
    - \usepackage{fullpage}
---

# Background

We show that aggregate usage (or, almost equivalently, the fraction using) and
the Gini coefficient of usage among users are both correlated with resistance.
To what degrees are these metrics comparable? Specifically, if we decrease one
marginal user by one claim (i.e., take a person with one claim to zero claims)
or a high-end user (e.g., take the largest user down by one claim), then how do
these two metrics change?

# Results

## Aggregate usage and fraction using

This are both simple. Say there are $B$ beneficiaries with $C$ total claims.
Then decreasing usage by one claim changes total usage by $1/C$ and changes the
fraction of nonzero users by $1/B$.

## Gini coefficient

Say there are $U$ beneficiaries with at least one claim, and $u_i$ is the usage
of the $i$-th user. Then the Gini coefficient of usage among users is
$$
G_+ \equiv \frac{\sum_i \sum_j \left| x_i - x_j \right|}{2 U C} = \frac{\sum_{i=0}^{U-1} \sum_{j=i+1}^U \left| x_i - x_j \right|}{U C}
$$
Write this as a fraction $\alpha/\beta$ so that
$$
\frac{\partial G_+}{\partial x_0} = \frac{\alpha'}{\beta} - \frac{\alpha \beta'}{\beta^2}.
$$
Consider the case where $x_0$ is the user with the most claims. Then
$$
\begin{aligned}
\alpha' &= \frac{\partial}{\partial x_0} \sum_{i=0}^{U-1} \sum_{j=i+1}^U \left| x_i - x_j \right| \\
  &= \frac{\partial}{\partial x_0} \sum_{j=1}^U \left( x_0 - x_j \right) \\
  &= U - 1
\end{aligned}
$$
and
$$
\beta' = \frac{\partial}{\partial x_0} UC = \frac{\partial}{\partial x_0} U \sum_i x_i = U
$$
so that
$$
\begin{aligned}
\frac{\partial G_+}{\partial x_0} &= \frac{U-1}{UC} - \frac{\alpha}{\beta} \frac{U}{UC} \\
  &= \frac{1 - 1/U}{C} - \frac{G_+}{U} \\
  &= \frac{1 - G_+ - 1/U}{C}.
\end{aligned}
$$
If $U \gg 1$, then this approaches $(1 - G_+)/C$, which is on the same order as
the change in total claims $1/C$ that came from changing the usage of a
marginal beneficiary.

## Discrete change

### Removing low consumer

Let the $x_i$ be ordered and let $x_1 = 1$. Then the question is what a new Gini coefficient looks like when removing a 1. Let $G_1$ be the "new" Gini coefficent and $G_0$ be the original:
$$
\begin{aligned}
G_1 + \frac{1}{(n-1)(S-1)} \sum_{j=2}^N (x_j - x_1) &= \frac{1}{(n-1)(S-1)} \sum_{i=1}^{N-1} \sum_{j=i+1}^N (x_j - x_i) \\
G_1 + \frac{S-n}{(n-1)(S-1)} &= \frac{nS}{(n-1)(S-1)} G_0 \\
G_1 - G_0 &= \frac{n + S - 1}{(n-1)(S-1)} G_0 - \frac{S-n}{(n-1)(S-1)}
\end{aligned}
$$
For $n \gg 1$, $S \gg 1$, we can throw out the 1's:
$$
G_1 - G_0 \approx \frac{S+n}{nS}G_0 - \frac{S-n}{nS}
$$

In our situation, many of our $x_i$ are ones, so $S$ is not very much larger
than $n$, and $G_0$ is not tiny, so the first term is larger than the first.
Applying the very crude estimate that $S \approx n$, then we get $\Delta \approx 2G_0/n$.

### Removing a top consumer

Let the $x_i$ be decreasing-ordered and let $x_1$ be the top donor. We ask what a new Gini looks like when change $x_1$ to $x_1' = x_1 - 1$.

# Discussion

This suggests an easier way to do the analysis: let $C/N$ be the total number
of claims per capita. We're interested in breaking down that $C/N$ into two
components: small-time users and big-time users. So write $C = A + B$, where
$A$ is the sum of all claims from users with a small number of claims (e.g., 0
or 1) and $B$ are the claims from all the big-time users (in this example, 2 or
greater).  Now $A/N$ and $B/N$ are on the same scale!
