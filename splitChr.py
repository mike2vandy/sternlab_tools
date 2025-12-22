#! /usr/bin/env python

import sys

with open(sys.argv[1]) as f: 
    for line in f:
        if 'chrUn' not in line:
            line = line.strip()
            fields = line.split()
            chrm, length = fields[0], int(fields[1])
            if length > 5000000:
                fifth = length / 5
                fifth = int(fifth)
                print(f"{chrm}:0-{fifth}")
                print(f"{chrm}:{fifth}-{fifth*2}")
                print(f"{chrm}:{fifth*2}-{fifth*3}")
                print(f"{chrm}:{fifth*3}-{fifth*4}")
                print(f"{chrm}:{fifth*4}-{length}")
            else:
                print(f"{chrm}:0-{length}")
