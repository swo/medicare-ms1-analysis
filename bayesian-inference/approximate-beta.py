#!/usr/bin/env python3

# How good is an approximation of the beta function?

import numpy as np
import scipy, scipy.integrate, scipy.optimize
from scipy.special import beta as beta_f

def est(a, b):
    approx = np.sqrt(2 * np.pi) * (a ** (a - 0.5) * b ** (b - 0.5)) / ((a + b) ** (a + b - 0.5))
    exact = beta_f(a, b)
    err = (exact - approx) / exact
    return exact, approx, err

vals = [1, 10, 100]
print('a', 'b', 'exact', 'approx', 'fractional error')
for a in vals:
    for b in [x for x in vals if x <= a]:
        print(a, b, est(a, b))
