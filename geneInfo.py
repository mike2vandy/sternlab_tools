#! /usr/bin/env python

import sys

print("geneName\tbiotype")
with open(sys.argv[1]) as f:
  for line in f:
    line = line.strip()
    if not line.startswith('#'):
      fields = line.split('\t')
      anno = fields[2]
      if anno == "gene":
        desc = fields[8].split(';')
        name = [i for i in desc if 'gene_id' in i][0].replace('"','').split()[1]
        biotype = [i for i in desc if 'gene_biotype' in i][0].replace('"','').split()[1]
        print(f"{name}\t{biotype}")
