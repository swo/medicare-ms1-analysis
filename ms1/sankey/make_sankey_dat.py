#!/usr/bin/env python3

years = [str(x) for x in range(2011, 2014+1)]

dat = {}
with open('tbl_cohorts.tsv') as f:
    header = next(f)
    for line in f:
        year_list, n = line.rstrip().split('\t')
        dat[year_list] = int(n)

n_bene = sum([dat[x] for x in dat])

with open('node.tsv', 'w') as f:
    def myprint(*x):
        print(*x, sep='\t', file=f)

    myprint('ID', 'x', 'y')
    for i, year in enumerate(years):
        myprint(year + 'in', i, 0)
        myprint(year + 'out', i, 1)

with open('edge.tsv', 'w') as f:
    def myprint(*x):
        print(*x, sep='\t', file=f)

    myprint('ID', 'N1', 'N2', 'Value')
    for year, next_year in zip(years, years[1:]):
        stay_in = sum([dat[x] for x in dat if year in x and next_year in x])
        leave = sum([dat[x] for x in dat if year in x and next_year not in x])
        enter = sum([dat[x] for x in dat if year not in x and next_year in x])
        stay_out = n_bene - (stay_in + leave + enter)

        myprint(year + 'stayin', year + 'in', next_year + 'in', stay_in)
        myprint(year + 'leave', year + 'in', next_year + 'out', leave)
        myprint(year + 'enter', year + 'out', next_year + 'in', enter)
        myprint(year + 'stayout', year + 'out', next_year + 'out', stay_out)
