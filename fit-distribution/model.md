---
title: 'Broken power law'
author: 'Scott Olesen'
header-includes:
    - \usepackage{fullpage}
---

The probability distribution function is
$$
f(x; \xi, \alpha_1, \alpha_2) = C(\xi, \alpha_1, \alpha_2) \times \begin{cases}
  x^{-\alpha_1} & \text{if } x < \xi \\
  \xi^{\alpha_2-\alpha_1} x^{-\alpha_2} & \text{if } x \geq \xi
\end{cases}
$$
Note that if $x = \xi$ then the two cases reduce to one another.

We find the normalization constant:
$$
\begin{aligned}
\frac{1}{C} &= \sum_{x=1}^\infty f(x) \\
 &= \sum_{x=1}^\xi f(x) + \sum_{x=\xi}^\infty f(x) \\
 &= \sum_{x=1}^\xi x^{-\alpha_1} + \xi^{\alpha_2 -\alpha_1} \sum_{x=\xi + 1}^\infty x^{-\alpha_2}
\end{aligned}
$$
Note that the second term is related to the zeta function:
$$
\zeta(k) \equiv \sum_{x=1}^\infty x^{-k}
$$
such that
$$
\sum_{x=\xi+1}^\infty x^{-\alpha_2} =
  \sum_{x=1}^\infty x^{-\alpha_2} - \sum_{x=1}^\xi x^{-\alpha_2} =
  \zeta(\alpha_2) - \sum_{x=1}^\xi x^{-\alpha_2}.
$$
Thus, the normalization constant is defined by:
$$
\frac{1}{C} = \sum_{x=1}^\xi x^{-\alpha_1} + \xi^{\alpha_2-\alpha_1} \times \left[
  \zeta(\alpha_2) - \sum_{x=1}^\xi x^{-\alpha_2}
\right]
$$

The likelihood of the data is $\mathcal{L}(X) = \prod_i f(x_i)$. Let's imagine that there are
$n_x$ data points for at each integer $x \geq 1$ so that
$$
\mathcal{L}(X) = C \times \prod_{x=1}^\xi \left[ x^{-\alpha_1} \right]^{n_i}
  \times \prod_{x=\xi+1}^\infty \left[ \xi^{\alpha_2-\alpha_1} x^{-\alpha_2} \right]^{n_i}.
$$
The log-likelihood is
$$
\mathcal{\ell}(X) = -\log \frac{1}{C} + \sum_{x=1}^\xi \left[ -n_x \alpha_1 \log x \right]
  + \sum_{x=\xi+1}^\infty \left[ n_x(\alpha_2-\alpha_1)\xi - n_x \alpha_2 \log x \right]
$$
