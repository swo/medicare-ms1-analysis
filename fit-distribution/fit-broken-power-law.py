#!/usr/bin/env python3

import numpy as np
import scipy.optimize, scipy.special

with open('n_claims.tsv') as f:
    header = next(f)
    rank_count = [list(map(int, line.rstrip().split())) for line in f]

def riemann_zeta(x):
    return 1 + scipy.special.zetac(x)

def inverse_broken_power_law_constant(rank_threshold, alpha1, alpha2):
    '''1/C, where C is the normalization constant'''
    term1 = sum([x ** (-alpha1) for x in range(1, rank_threshold + 1)])
    xaa = rank_threshold ** (alpha2 - alpha1)
    term2 = riemann_zeta(alpha2) - sum([x ** (-alpha2) for x in range(1, rank_threshold + 1)])
    result = term1 + xaa * term2

    if result <= 0:
        print('t1',term1)
        print('xaa',xaa)
        print('t2',term2)
        raise ValueError('bad normalization constant {} for rt={} a1={} a2={}'.format(result, rank_threshold, alpha1, alpha2))
    else:
        return result

def broken_power_law_density(rank_threshold, alpha1, alpha2, x):
    c = 1.0 / inverse_broken_power_law_constant(rank_threshold, alpha1, alpha2)

    if alpha1 <= 1.0 or alpha2 <= 1.0:
        raise ValueError('bad alphas: {}, {}'.format(alpha1, alpha2))

    if x < 1:
        raise ValueError
    elif x < rank_threshold:
        result = c * x ** (-alpha1)
    else:
        result = c * (rank_threshold ** (alpha2 - alpha1)) * (x ** (-alpha2))

    if result <= 0:
        print('c', c)
        raise ValueError('bad density {} for rt={} a1={} a2={} x={}'.format(result, rank_threshold, alpha1, alpha2, x))
    else:
        return result

# def log_likelihood(rank_count, rank_threshold, alpha1, alpha2):
#     term1 = -np.log(inverse_broken_power_law_constant(rank_threshold, alpha1, alpha2))

#     low_rank_count = [x for x in rank_count if x[0] < rank_threshold]
#     high_rank_count = [x for x in rank_count if x[0] >= rank_threshold]
#     term2 = sum([-n * alpha1 * np.log(x) for x, n in low_rank_count])
#     term3 = sum([n * (alpha2 - alpha1) * rank_threshold - n * alpha2 * np.log(x) for x, n in high_rank_count])

#     return term1 + term2 + term3

def log_likelihood(rank_count, rank_threshold, alpha1, alpha2):
    return sum([n * np.log(broken_power_law_density(rank_threshold, alpha1, alpha2, x)) for x, n in rank_count])

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

print(broken_power_law_density(30, 1.6369327492532373, 17.11981591110122, 50))
print(broken_power_law_density(30, 1.6369327492532373, 17.11981591110122, 1))

xis = [15, 20, 25, 30]
results = {xi: broken_at(xi) for xi in xis}
best_xi = max(xis, key=lambda xi: results[xi].fun)
print('best xi', best_xi)
print(results[best_xi])
