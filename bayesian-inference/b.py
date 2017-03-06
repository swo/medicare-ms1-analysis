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


def logistic(par_l, par_x0, par_k, x):
    return par_l / (1 + np.exp(-par_k * (x - par_x0)))

def log_likelihood(par, states):
    par_l, par_x0, par_k = par
    mean_resistances = [logistic(par_l, par_x0, par_k, s.consumption) for s in states]
    assert all([m > 0 for m in mean_resistances])

    def term_ij(s, n, m, strength):
        a = m * strength
        b = (1 - m) * strength
        assert m > 0 and m < 1
        assert a > 0
        assert b > 0
        return beta_f(a + s, b + n - s) / beta_f(a, b)

    def term_i(mean, state):
        f = lambda s: np.prod([term_ij(h.n_resistant, h.n_isolates, mean, s) for h in state.hospitals])
        y, err = scipy.integrate.quad(f, 0, np.inf)

        if y == 0:
            print('mean', mean)
            print('nr, ni', [(h.n_resistant, h.n_isolates) for h in state.hospitals])

        assert y > 0
        return y

    return sum([np.log(term_i(m, s)) for m, s in zip(mean_resistances, states)])

def mle(data):
    f = lambda x: -log_likelihood(x, data)
    x0 = np.array([0.5, 0.5, 0.5])

    res = scipy.optimize.minimize(f, x0, options={'disp': True})
    return res

n_states = 10

par_l = 0.25
par_x0 = 0.0
par_k = 1.0

states = []
for i in range(n_states):
    consumption = np.random.uniform(0, 100)
    mean = logistic(par_l, par_x0, par_k, consumption)

    strength = np.random.uniform(2, 20)
    var = mean * (1 - mean) / (1 + strength)
    alpha = mean * strength
    beta = (1 - mean) * strength

    n_hospitals = np.random.random_integers(low=1, high=20)
    hospitals = []
    for i in range(n_hospitals):
        resistance = np.random.beta(alpha, beta)
        n_isolates = np.random.random_integers(1, 500)
        n_resistant = np.random.binomial(n_isolates, resistance)
        hospitals.append(Hospital(n_isolates, n_resistant))

    states.append(State(consumption, hospitals))

res = mle(states)

print(res)
