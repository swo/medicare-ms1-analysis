#!/usr/bin/env python3

import numpy as np
import scipy, scipy.optimize

dat = {}
with open('ecoli_quin.tsv') as f:
    header = next(f)
    for line in f:
        state, nr, ni = line.rstrip().split('\t')

        if ni == 'NA':
            continue

        nr = int(float(nr))
        ni = int(ni)
        if state not in dat:
            dat[state] = []

        dat[state].append((nr, ni))

def log_likelihood(mu, nu, data):
    alpha = mu * nu
    beta = (1.0 - mu) * nu
    return sum([scipy.special.betaln(alpha + s, beta + n - s) for s, n in data]) - \
            len(data) * scipy.special.betaln(alpha, beta)

def mle(data):
    prior = lambda m, n: np.exp(-n)
    f = lambda x: -prior(x[0], x[1]) * log_likelihood(x[0], x[1], data)

    total_r = sum([x[0] for x in data])
    total_i = sum([x[1] for x in data])
    grand_mean = total_r / total_i

    x0 = np.array([grand_mean, 1.0])

    res = scipy.optimize.minimize(f, x0, bounds=[(0, 1), (0, np.inf)], options={'disp': False})
    return res

results = {}
for state in dat:
    res = mle(dat[state])
    mu, nu = res.x
    results[state] = (mu, nu)

print('state', 'mu', 'nu', sep='\t')
for state in sorted(results.keys()):
    print(state, results[state][0], results[state][1], sep='\t')
