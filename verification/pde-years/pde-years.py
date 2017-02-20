#!/usr/bin/env python3

import re, argparse, sys

def parse_line(line):
    date = line.rstrip().split('\t')[1]
    assert re.match('\d{2}[A-Z]{3}\d{2}', date)
    day = date[0:2]
    month = date[2:5]
    year = date[5:9]
    return (year, month, day)

def summarize(lines):
    data = {}

    for line_i, line in enumerate(lines):
        if line_i == 0:
            # this is the header
            continue
        
        date = parse_line(line)
        if date not in data:
            data[date] = 0

        data[date] += 1

    return data

if __name__ == '__main__':
    p = argparse.ArgumentParser(description='parse PDE files and summarize dates')
    p.add_argument('input', type=argparse.FileType('r'), nargs='?', default=sys.stdin)
    args = p.parse_args()

    result = summarize(args.input)
    print('year', 'month', 'day', 'count', sep='\t')
    for k in sorted(result.keys()):
        print(k[0], k[1], k[2], result[k], sep='\t')
