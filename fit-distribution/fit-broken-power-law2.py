#!/usr/bin/env python3

import numpy as np
import scipy.optimize, scipy.special, scipy.stats

with open('n_claims.tsv') as f:
    header = next(f)
    rank_count = [list(map(int, line.rstrip().split())) for line in f]

def log_likelihood(rank_count, alpha):
    return sum([n * scipy.stats.zipf.logpmf(x, alpha) for x, n in rank_count])

def protect_alpha(alpha):
    return 1.0 + alpha ** 2

def optimize_power_law(rank_count):
    x0 = 2.0
    f = lambda alpha: -log_likelihood(rank_count, protect_alpha(alpha))
    res = scipy.optimize.minimize(f, x0)
    print('optimal alpha', protect_alpha(res.x))
    return res

def objective_function(rank_threshold, x):
    alpha1, alpha2 = x
    alpha1 = 1.0 + alpha1 ** 2
    alpha2 = 1.0 + alpha2 ** 2
    return -log_likelihood(rank_count, rank_threshold, alpha1, alpha2)

def broken_at(rank_threshold):
    alpha1_0 = 1.0
    alpha2_0 = 1.0
    x0 = np.array([alpha1_0, alpha2_0])
    f = lambda x: objective_function(rank_threshold, x)
    res = scipy.optimize.minimize(f, x0)
    return res

#print(broken_power_law_density(30, 1.6369327492532373, 17.11981591110122, 50))
#print(broken_power_law_density(30, 1.6369327492532373, 17.11981591110122, 1))

#xis = [15, 20, 25, 30]
#results = {xi: broken_at(xi) for xi in xis}
#best_xi = max(xis, key=lambda xi: results[xi].fun)
#print('best xi', best_xi)
#print(results[best_xi])

rc = [x for x in rank_count if x[0] < 11]
print(optimize_power_law(rc))
rc = [x for x in rank_count if x[0] >= 11]
print(optimize_power_law(rc))
