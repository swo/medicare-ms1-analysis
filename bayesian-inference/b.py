#!/usr/bin/env python3

import numpy as np
import scipy, scipy.integrate, scipy.optimize
from scipy.special import beta as beta_f

class State:
    def __init__(self, consumption, hospitals):
        self.consumption = consumption
        self.hospitals = hospitals


class Hospital:
    def __init__(self, n_isolates, n_resistant):
        self.n_isolates = n_isolates
        self.n_resistant = n_resistant


def alpha_coeff(m, v):
    x = (m ** 2 - m ** 3 - m * v) / v
    print(m, v)
    assert x > 0
    return x

def beta_coeff(m, v):
    x = (m - 1) * (m ** 2 - m + v) / v
    print(m, v)
    assert x > 0
    return x

def log_likelihood(slope, states):
    mean_resistances = [s.consumption * slope for s in states]
    print(mean_resistances)

    def term_ij(s, n, m, v):
        a = alpha_coeff(m, v)
        b = beta_coeff(m, v)
        assert m > 0
        assert v > 0
        return beta_f(a + s, b + n - s) / beta_f(a, b)

    def term_i(mean, state):
        f = lambda v: np.prod([term_ij(h.n_resistant, h.n_isolates, mean, v) for h in state.hospitals])
        #f = lambda v: np.prod([term_ij(s, n, mean, v) for s, n in counts])
        y, err = scipy.integrate.quad(f, 0, mean * (1 - mean))
        return y

    return sum([np.log(term_i(m, s)) for m, s in zip(mean_resistances, states)])

def mle(data):
    f = lambda x: -log_likelihood(x, data)
    x0 = 1.0

    res = scipy.optimize.minimize(f, x0, bounds=[(1e-6, 20)], options={'disp': True})
    return res

n_states = 10
slope = 0.25

data = []
for i in range(n_states):
    mean = np.random.uniform(0, 1)
    consumption = mean / slope
    strength = np.random.uniform(2, 20)
    var = mean * (1 - mean) / (1 + strength)
    alpha = alpha_coeff(mean, var)
    beta = beta_coeff(mean, var)

    n_hospitals = np.random.random_integers(low=1, high=20)
    hospitals = []
    for i in range(n_hospitals):
        resistance = np.random.beta(alpha, beta)
        n_isolates = np.random.random_integers(1, 500)
        n_resistant = np.random.binomial(n_isolates, resistance)
        hospitals.append(Hospital(n_isolates, n_resistant))

    data.append(State(consumption, hospitals))

print(data)

# state_means = np.random.uniform(0, 1, size=n_states)
# state_strengths = np.random.uniform(2, 20, size=n_states)
# state_variances = [m * (1 - m) / (1 + n) for m, n in zip(state_means, state_strengths)]
# consumptions = [m / slope for m in state_means]
# n_hospitals = np.random.random_integers(low=1, high=20, size=n_states)
# n_isolates = [np.random.random_integers(low=1, high=500, size=nh) \
#         for nh in n_hospitals]
# hospital_resistances = [np.random.beta(alpha_coeff(m, v), beta_coeff(m, v), size=nh) \
#         for m, v, nh in zip(state_means, state_variances, n_hospitals)]
# n_resistant = [[np.random.binomial(ni, r) for ni, r in zip(nis, rs)] \
#         for nis, rs in zip(n_isolates, hospital_resistances)]


# counts = list(zip(n_resistant, n_isolates))
# print(counts)
# res = mle(consumptions, counts)

res = mle(data)

print(res)

# print('mle ab', res.x)
# print('mle mean', res.x[0] / sum(res.x))
# print('mle var', bvar(*res.x))
# print('true ab', alpha, beta)
# print('true mean', alpha / (alpha + beta))
# print('true var', bvar(alpha, beta))
# print('ns', ns)
# print('ps', ps)
# print('data', data)
