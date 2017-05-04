#!/usr/bin/env python3

# Compute power of Spearman's rho test for given size and rho
# Draw pairs from a bivariate normal with that rho. Then compute
# the number of cases below p=0.05.

import scipy.stats

def gini(f, mean, n):
    return 1.0 / mean * sum([f(i) * f(j) * (i - j) for i in range(1, n) for j in range(0, i)])
    #return sum([f(i) * f(j) * abs(i - j) for i in range(0, n) for j in range(0, n)]) / (2 * mean)

def gini_eps(f, mean, eps=1e-12, verbose=False):
    total = 0.0
    term = 0.0
    i = 1
    while i < 10 or i < mean or term > eps:
        if verbose:
            print('term i', i, term, 'total', total)

        term = 1.0 / mean * f(i) * sum([f(j) * (i - j) for j in range(0, i - 1)])
        total += term
        i += 1

    return total

def gini_nbinom(r, p, eps=1e-12):
    f = lambda k: scipy.stats.nbinom.pmf(k, r, p)
    mu = p * r / (1 - p)
    # print('mu', mu)
    return gini_eps(f, mu, eps=eps)
    #return gini(f, mu, 40)

print('r', 'p', 'gini', 'nzgini', sep='\t')
# for r in [0.1, 0.5, 1, 2, 5]:
#     for p in [0.1, 0.25, 0.5, 0.75, 0.99]:
for r in [0.25, 0.5, 0.75]:
    for p in [0.275, 0.3, 0.325]:
        f = lambda k: scipy.stats.nbinom.pmf(k, r, p)
        mu = scipy.stats.nbinom.stats(r, p, moments='m')

        def f2(k):
            if k == 0:
                return 0
            else:
                return f(k) / (1.0 - f(0))

        mu2 = mu / (1.0 - f(0))
        #print(r, p, gini_nbinom(r, p), sep='\t')
        print(r, p, gini(f, mu, 100), gini(f2, mu2, 100), sep='\t')

# r = 0.1
# p = 0.1
# f = lambda k: scipy.stats.nbinom.pmf(k, r, p)
# mu = scipy.stats.nbinom.stats(r, p, moments='m')
# print('gini', gini(f, mu, 100))
#print([f(i) for i in range(10)])

#mu = 3.0
#f = lambda k: scipy.stats.poisson.pmf(k, mu)
#print(gini(f, mu, 10))

# print(gini_nbinom(1, 0.1))
# print(scipy.stats.nbinom.pmf(0, 1, 0.1))
