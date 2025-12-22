#! /usr/bin/env python

from cyvcf2 import VCF
import sys

iupac = {"AG": "R", "CT":"Y", "AC":"M", "GT":"K",
         "CG": "S", "AT": "W"}


vcf_file = VCF(sys.argv[1])

samples = vcf_file.samples

print(samples)
seqs = {i: '' for i in samples}

for record in vcf_file:
  ref, alt = record.REF, record.ALT[0]
  if '*' not in ref and '*' not in alt:
    var_geno = dict(zip(samples, record.genotypes))
    for sample, gt in var_geno.items():
      a1, a2 = gt[0], gt[1]
      if a1 == 0 and a2 == 0:
        seqs[sample] += ref
      elif a1 == 1 and a2 == 1:
        seqs[sample] += alt
      else:
        code = iupac[''.join(sorted([ref,alt]))]
        seqs[sample] += code

for i, j in seqs.items():
  print(f">{i}\n{j}")
 
