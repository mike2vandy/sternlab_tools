#! /usr/bin/env python

import sys
import pandas as pd
from cyvcf2 import VCF
import openpyxl


vcf_file = VCF(sys.argv[1])

meta = pd.read_csv(sys.argv[2])

region = sys.argv[3]

out = sys.argv[4]

samples = vcf_file.samples

trans = {'0/0':0, '0/1':1, '1/1':2}

rows = []
for record in vcf_file(region):
  genos = dict(zip(samples, record.genotypes))
  for sample, gt in genos.items():
    new_gt = '/'.join(map(str,gt[:2]))
    if '-' in new_gt:
      best = -1
    else:
      best = trans[new_gt]
    
    rows.append({"sample" : sample, "n_alleles" : best}) 
   
df = pd.DataFrame(rows)

main = pd.merge(df, meta, left_on = "sample", right_on = "DogID", how = 'inner')

main_sub = main[["sample", "n_alleles", "Breed"]]

main_sub.to_excel(sys.argv[5])

main_sub = main_sub[main_sub['n_alleles'] != -1]

breed_af = main_sub.groupby("Breed").agg(n_dogs = ('n_alleles', 'count'), total_alleles = ('n_alleles', 'sum'))

breed_af['AF'] = round(breed_af['total_alleles'] / (breed_af['n_dogs']*2), 3)

breed_af.to_excel(out)
