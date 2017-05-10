#!/usr/bin/env python3

# Compute p-values for Spearman correlation using the Fisher transformation
# as explained on wikipedia.

import numpy as np
import scipy.stats

def max_sdi2(n):
    '''Maximum sum of d_i^2 (i.e., for which r=0)'''
    return n * (n ** 2 - 1) // 6

def r(sdi2, n):
    '''Correlation coefficient'''
    return 1.0 - 6 * sdi2 / (n * (n ** 2 - 1))

def z_score(sdi2, n):
    '''Approximate z score'''
    if sdi2 == 0:
        return np.inf
    else:
        return np.sqrt((n - 3) / 1.06) * np.arctanh(r(sdi2, n))

def p_value(sdi2, n):
    '''Approximate p value'''
    z = z_score(sdi2, n)
    return scipy.stats.norm.sf(z) + scipy.stats.norm.cdf(-z)

print('n', 'sdi2', 'r', 'z', 'p', sep='\t')
for n in [19, 30, 35, 40]:
    for sdi2 in range(0, max_sdi2(n) + 1, 100):
        print(n, sdi2, r(sdi2, n), z_score(sdi2, n), p_value(sdi2, n), sep='\t')
