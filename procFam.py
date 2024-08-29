#! /usr/bin/env python

import sys, vcf

#best for dog snv vcf from WAGS
def reduceCSQ(csqs):
    #columns: variant,consequence,impact,gene,transcriptID,biotype,exon,hgvsc,hgvsp,cDNA_pos,cds_pos,aa_pos,aas,codons,variant_class,proteinID
    mainCols = [0, 1, 2, 4, 5, 6, 7, 8, 9, 12, 13, 14, 15, 16, 19, 21, 30]
    updated = []
    for record in csqs:
        record = record.split('|')
        result = ','.join([record[i] for i in mainCols])
        updated.append(result)

    return(updated)

sampleNames = ['5370_Lincoln', '5387_Whiskey', '5377_Tenney', '5378_Natasha', '5380_BB', '5384_Barley', '5388_Grizzly']

lincoln = '5370_Lincoln'
father = '5387_Whiskey'
siblings = ['5377_Tenney',
    '5378_Natasha',
    '5380_BB',
    '5384_Barley',
    '5388_Grizzly']

vcfR = vcf.Reader(filename = sys.argv[1])

header = (
  "chrm,pos,ref,alt,filter,AF,allele,"
  "Consequence,IMPACT,Gene,Feature_type,Feature,BIOTYPE,EXON,INTRON,cDNA_position,"
  "CDS_position,Protein_position,Amino_acids,Codons,STRAND,VARIANT_CLASS,protein_id,"
  f"{','.join(sampleNames)}")

print(header)
for record in vcfR:
  chrm, pos, ref, alt = record.CHROM, record.POS, record.REF, '/'.join(map(str, record.ALT))
  if record.FILTER:
    filt = record.FILTER[0]
  else:
    filt = 'PASS'
  info = record.INFO
  AF = info.get('AF', ['NA'])[0]
  main = ''
  dad = ''
  sibs = []
  for i in record.samples: 
    sam = i.sample
    gt = i['GT'].replace('|', '/')
    
    if sam == lincoln and gt == '1/1':
      main = gt
    elif sam == father and gt == '0/1':
      dad = gt
    elif sam in siblings and gt != '1/1':
      sibs.append(gt)
    else:
      pass
  if main and dad and len(sibs) == 5:
    csq = reduceCSQ(info.get('CSQ', [])) 
    for consequence in csq:
      print(f"{chrm},{pos},{ref},{alt},{filt},{AF},{consequence},{main},{dad},{','.join(sibs)}") 
