#!/usr/bin/env python3

import fileinput

class DxReader:
    top_abx = [None, 'azithromycin', 'ciprofloxacin', 'amoxicillin', 'cephalexin',
        'trimethoprim/sulfamethoxazole', 'levofloxacin', 'amoxicillin/clavulanate',
        'doxycycline', 'nitrofurantoin', 'clindamycin']

    out_abx = ['.n', 'azithromycin', 'ciprofloxacin', 'amoxicillin', 'cephalexin',
        'trimethoprim/sulfamethoxazole', 'levofloxacin', 'amoxicillin/clavulanate',
        'doxycycline', 'nitrofurantoin', 'clindamycin']

    def __init__(self):
        self.bene = None
        self.dat = {}
        self.first_line = True

        print('bene_id', 'dx_cat', 'n_encounters', *self.out_abx[1:], sep='\t')

    def parse_line(self, line):
        fields = line.rstrip().split('\t')
        bene, enc, dx = fields[0:3]

        if len(fields) == 4:
            abx = fields[3]
        else:
            abx = None

        return (bene, enc, dx, abx)

    def process_line(self, line):
        bene, enc, dx, abx = self.parse_line(line)

        if self.first_line:
            self.bene = bene
            self.dat = {}
            self.add_to_dat(enc, dx, abx)
            self.first_line = False
        else:
            if bene != self.bene:
                self.report()
                self.dat = {}
                self.bene = bene

            self.add_to_dat(enc, dx, abx)

    def add_to_dat(self, enc, dx, abx):
        if enc in self.dat:
            self.dat[enc]['.dx'].append(dx)

            if abx in self.dat[enc]:
                self.dat[enc][abx] += 1
            else:
                self.dat[enc][abx] = 1
        else:
            self.dat[enc] = {'.dx': [dx], abx: 1}

    def report(self):
        dxs = set()
        for enc in self.dat:
            dxs.update(self.dat[enc]['.dx'])

        out = {dx: {'.n': 0} for dx in dxs}

        for enc in self.dat:
            for dx in self.dat[enc]['.dx']:
                out[dx]['.n'] += 1
                for abx in self.dat[enc].keys():
                    if abx == '.dx':
                        continue

                    if abx not in self.top_abx:
                        abx = 'other'

                    if abx is not None:
                        if abx in out[dx]:
                            out[dx][abx] += 1
                        else:
                            out[dx][abx] = 1

        for dx in out:
            out_fields = [self.bene, dx] + [out[dx][abx] if abx in out[dx] else 0 for abx in self.out_abx]
            print(*out_fields, sep='\t')

if __name__ == '__main__':
    r = DxReader()
    for line in fileinput.input():
        r.process_line(line)

    r.report()
