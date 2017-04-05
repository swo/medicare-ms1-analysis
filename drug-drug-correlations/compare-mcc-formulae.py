#!/usr/bin/env python3

'''
Confirm that the two ways to compute the MCC that I read were equivalent.
'''

import numpy as np
import scipy.stats

def naive(true, pred):
    tp = len([1 for t, p in zip(true, pred) if t and p])
    tn = len([1 for t, p in zip(true, pred) if not t and not p])
    fp = len([1 for t, p in zip(true, pred) if not t and p])
    fn = len([1 for t, p in zip(true, pred) if t and not p])

    denom = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)

    if denom == 0:
        return 0
    else:
        return (tp * tn - fp * fn) / np.sqrt(denom)

def scikit(true, pred):
    mean_true = np.average(true)
    mean_pred = np.average(pred)
    u_true = np.array(true) - mean_true
    u_pred = np.array(pred) - mean_pred

    cov = np.average(u_true * u_pred)
    var_true = np.average(u_true ** 2)
    var_pred = np.average(u_pred ** 2)

    denom = var_true * var_pred

    if denom == 0:
        return 0
    else:
        return cov / np.sqrt(var_true * var_pred)

print(naive([0,1,0,1,0], [1,1,1,0,0]))
assert False

for i in range(int(1e3)):
    true = scipy.stats.binom.rvs(n=1, p=0.5, size=10)
    pred = scipy.stats.binom.rvs(n=1, p=0.5, size=10)

    x = naive(true, pred)
    y = scikit(true, pred)

    if abs(x - y) > 1e-12:
        raise RuntimeError
    else:
        print(x)
