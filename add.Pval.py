#! /usr/bin/env python

import sys

pvals = {}
with open(sys.argv[1]) as f:
  for line in f:
    fields = line.strip().split()
    site, pval = fields[1], fields[10]
    pvals[site] = pval

with open(sys.argv[2]) as f:
  for line in f:
    line = line.strip()
    if line.startswith('chrm'):
      print(f"{line},Wald_p")
    else:
      fields = line.split(',')
      index = f"{fields[0]}_{fields[1]}"
      if index in pvals:
        print(f"{line},{pvals[index]}")
      else:
        print(f"{line},NA")
