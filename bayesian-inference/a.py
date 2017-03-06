#!/usr/bin/env python3

import numpy as np
import scipy, scipy.optimize

def bvar(a, b):
    return a * b / ((a + b) ** 2 * (a + b + 1))

def log_likelihood(alpha, beta, data):
    return sum([scipy.special.betaln(si + alpha, ni - si + beta) for si, ni in data]) - \
            len(data) * scipy.special.betaln(alpha, beta)

def mle(data):
    f = lambda x: -log_likelihood(x[0], x[1], data)
    x0 = np.array([1.0, 1.0])

    res = scipy.optimize.minimize(f, x0, bounds=[(1e-6, np.inf)] * 2, options={'disp': False})
    return res

ns = np.random.random_integers(low=1, high=500, size=10)

print('true_alpha', 'true_beta', 'true_mean', 'true_var', 'mle_alpha', 'mle_beta', 'mle_mean', 'mle_var', sep='\t')
for i in range(100):
    alpha, beta = np.random.uniform(low=1e-6, high=10, size=2)
    ps = np.random.beta(alpha, beta, size=len(ns))
    data = [(np.random.binomial(n, p), n) for n, p in zip(ns, ps)]
    res = mle(data)

    print(alpha, beta, alpha / (alpha + beta), bvar(alpha, beta), \
            res.x[0], res.x[1], res.x[0] / sum(res.x), bvar(*res.x), sep='\t')

# print('mle ab', res.x)
# print('mle mean', res.x[0] / sum(res.x))
# print('mle var', bvar(*res.x))
# print('true ab', alpha, beta)
# print('true mean', alpha / (alpha + beta))
# print('true var', bvar(alpha, beta))
# print('ns', ns)
# print('ps', ps)
# print('data', data)
