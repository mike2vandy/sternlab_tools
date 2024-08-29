#! /usr/bin/env python3

import sys
from collections import defaultdict

genesBound = defaultdict()
geneEx = defaultdict()

with open(sys.argv[1]) as f:
  for line in f:
    line = line.strip()
    fields = line.split('\t')
    chrm, start, stop = fields[0], int(fields[3]), int(fields[4])
    orient = fields[6]
    info = fields[8]
    gene = [i for i in info.split(';') if 'gene_id' in i][0].split()[1].replace('"', '')
    if fields[2] == 'transcript':
      genesBound[gene] = [chrm, start, stop, orient] 
    if fields[2] == 'exon':
      if gene in geneEx:
        geneEx[gene].append([start, stop])
      else:
        geneEx[gene] = [[start, stop]]

cdsOrder = {}
for i, j in geneEx.items():
  cdsOrder[i] = [[],[]]
  cdsCt = 1
  for k in j:
     genomeCt = k[0]
     while genomeCt <= k[1]:
       cdsOrder[i][0].append(genomeCt)
       cdsOrder[i][1].append(cdsCt)
       cdsCt += 1 
       genomeCt += 1

print("index,transPos,gene")
for i, j in cdsOrder.items():
  chrm, start, end, orient = genesBound[i]
  if genesBound[i][3] == '-':
    j[0].reverse()
    for k, l in zip(j[0], j[1]):
      print(f"{chrm}_{k},{l},{i}")
  else:
    for k, l in zip(j[0], j[1]):
      print(f"{chrm}_{k},{l},{i}")
   
''' 
    genePos = geneStart
    buff = 0
    while genePos > j[0][0]:
      buff += 1
      print(f"{chrm},{genePos},{buff},5-prime',{i}")#"i, genePos, buff, '5-prime')
      genePos -= 1
 
    cdsEnd = 0
    for k, l in zip(j[0], j[1]):
      print(f"{chrm},{k},{l+buff},cds,{i}")#i, k, l+buff, 'cds')
      cdsEnd = l+buff
      genePos = k

    while genePos > geneEnd:
      cdsEnd += 1
      genePos -= 1
      print(f"{chrm},{genePos},{cdsEnd},3-prime,{i}")#i, genePos, cdsEnd, geneEnd, '3-prime')

  else:
    geneStart, geneEnd = start, end
    genePos = geneStart

    buff = 0
    while genePos < j[0][0]:
      buff += 1
      print(f"{chrm},{genePos},{buff},5-prime,{i}")#i, genePos, buff, j[0][0], '5-prime')
      genePos += 1
    
    cdsEnd = 0
    for k, l in zip(j[0], j[1]):
      print(f"{chrm},{k},{l+buff},cds,{i}")# i, k, l+buff, 'cds')
      cdsEnd = l+buff 
      genePos = k 

    while genePos < geneEnd:
      cdsEnd += 1
      genePos += 1
      print(f"{chrm},{genePos},{cdsEnd},3-prime,{i}")#i, genePos, cdsEnd, geneEnd, '3-prime')
      
'''
