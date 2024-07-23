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

vcfR = vcf.Reader(filename = sys.argv[1])

sampleNames = vcfR.samples

#print()

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
  dogs = []
  for i in record.samples: 
    sam = i.sample
    gt = i['GT']
    dogs.append(gt)
  csq = reduceCSQ(info.get('CSQ', [])) 
  for consequence in csq:
    print(f"{chrm},{pos},{ref},{alt},{filt},{AF},{consequence},{','.join(dogs)}") 
