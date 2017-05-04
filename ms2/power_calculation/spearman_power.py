#!/usr/bin/env python3

# Compute power of Spearman's rho test for given size and rho
# Draw pairs from a bivariate normal with that rho. Then compute
# the number of cases below p=0.05.

import numpy as np
import scipy.stats

# Shouldn't matter what the exact sd's are
sd1 = 1.0
sd2 = 1.0

def power(size, rho, n_trials=1000):
    n_sig = 0
    for i in range(n_trials):
        cov = [[sd1 ** 2, rho * sd1 * sd2], [rho * sd1 * sd2, sd2 ** 2]]
        x = np.random.multivariate_normal([0, 0], cov, size=size)
        est_rho, p = scipy.stats.spearmanr(x)
        if p < 0.05:
            n_sig += 1

    return n_sig / n_trials

print('size', 'rho', 'power', sep='\t')
for size in [25, 30, 35, 40, 45]:
    for rho in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
        print(size, rho, power(size, rho), sep='\t')
