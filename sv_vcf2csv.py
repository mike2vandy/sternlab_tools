#! /usr/bin/env python

import sys, vcf

vcfR = vcf.Reader(filename = sys.argv[1])

sampleNames = vcfR.samples

#print()

print(f"chrm,pos,svTyp,length,AF,impact,consq,gene,trans,exon,biotyp,{','.join(sampleNames)}")
for record in vcfR:
  chrm, pos = record.CHROM, record.POS
  info = record.INFO
  af = info['AF'][0]
  length = info['SVLEN']
  svTyp = info['SVTYPE']
  dogs = []
  for i in record.samples: 
    sam = i.sample
    gt = i['GT'].replace('|','/')
    dogs.append(gt)
  try:
    cons = info['CSQ']
    for i in cons:
      tmp = i.split('|')
      impact, consq, gene, trans, biotyp, exon = tmp[1], tmp[2], tmp[4], tmp[6], tmp[7], tmp[9]
      print(f"{chrm},{pos},{svTyp},{length},{af},{impact},{consq},{gene},{trans},{exon},{biotyp},{','.join(dogs)}")
  except:
    pass
