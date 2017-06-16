#!/usr/bin/env python3

'''
Computations of the confidence intervals for relative risks using a normal
approximation for its logarithm. As per p.266ff in Practical Statistics for
Medical Research, D. G. Altman, Chapman and Hall, London, 1991
'''

import numpy as np, scipy.stats

def rrci(a, b, c, d, level):
    alpha = (1.0 - level) / 2
    var = scipy.stats.norm.isf(alpha)

    rr = (a / (a + c)) / (b / (b + d))
    log_rr = np.log(rr)
    se_log_rr = np.sqrt(1.0 / a - 1.0 / (a + c) + 1.0 / b - 1.0 / (b + d))

    deviation = var * se_log_rr

    return (np.exp(log_rr - deviation), rr, np.exp(log_rr + deviation))

def rrci_total(a, ac, b, bd, level):
    '''
    a : exposed & outcome
    ac : exposed
    b : not exposed & outcome
    bd: not exposed
    '''

    c = ac - a
    d = bd - b
    return rrci(a, b, c, d, level)

# Verification from p. 267 example
# print(rrci(2, 33, 14, 58, 0.9))

# Reuse the totals for the UTI drugs
n_women_uti = 237536
n_women_no_uti = 1318915
n_men_uti = 63137
n_men_no_uti = 858710

def display(drug, women_a, women_b, men_a, men_b, level=0.95):
    print('women', drug, rrci_total(women_a, n_women_uti, women_b, n_women_no_uti, level))
    print('men', drug, rrci_total(men_a, n_men_uti, men_b, n_men_no_uti, level))

display('cipro', 112984, 136368, 29263, 79416)
print()
display('tmp/smx', 59577, 81586, 15835, 47395)
print()
display('nitro', 65030, 38636, 8738, 5634)
print()

# Azithromycin uses respiratory complaints as denominators
print('women', 'azithro', rrci_total(177843, 402523, 109681, 1153928, 0.95))
print('men', 'azithro', rrci_total(84145, 208697, 54134, 713150, 0.95))

# Cipro, with the new data
print('women', 'cipro', 'newdata', rrci(122298, 127054, 138795, 1168304, 0.95))
