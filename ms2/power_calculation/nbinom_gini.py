#!/usr/bin/env python3

# Compute power of Spearman's rho test for given size and rho
# Draw pairs from a bivariate normal with that rho. Then compute
# the number of cases below p=0.05.

import scipy.stats

def gini(f, mean, n):
    return 1.0 / mean * sum([f(i) * f(j) * (i - j) for i in range(1, n) for j in range(0, i - 1)])

def gini_eps(f, mean, eps=1e-12, verbose=False):
    total = 0.0
    term = 0.0
    i = 1
    while i < mean or term > eps:
        if verbose:
            print('term i', i, term, 'total', total)

        term = 1.0 / mean * f(i) * sum([f(j) * (i - j) for j in range(0, i - 1)])
        total += term
        i += 1

    return total

def gini_nbinom(p, mu):
    r = (1 - p) * mu / p
    print('p', p, 'mu', mu, 'r', r)
    f = lambda k: scipy.stats.nbinom.pmf(k, r, p)
    n = 40
    return (gini(f, mu, n), gini(f, mu, n + 1))

def gini_nbinom_eps(p, mu, eps=1e-12):
    r = (1 - p) * mu / p
    print('p', p, 'mu', mu, 'r', r)
    f = lambda k: scipy.stats.nbinom.pmf(k, r, p)
    return gini_eps(f, mu, eps=eps)

print(gini_nbinom(0.5, 10))
print(gini_nbinom_eps(0.5, 10))

