#!/usr/bin/env python

import os
import sys
from collections import Counter
#from collections import defaultdict

import vcf

def reduceCSQ(csqs):
    #columns: variant,consequence,impact,gene,transcriptID,biotype,exon,hgvsc,hgvsp,cDNA_pos,cds_pos,aa_pos,aas,codons,variant_class,proteinID
    mainCols = [0, 1, 2, 4, 5, 6, 7, 8, 9, 12, 13, 14, 15, 16, 19, 21, 30]
    updated = []
    for record in csqs:
        record = record.split('|')
        result = ','.join([record[i] for i in mainCols])
        updated.append(result)
    return(updated)

def calcAF(groupList, minor_alt):
    #calculate AF
    alleles = [allele for genotype in groupList for allele in genotype.split('/') if allele != '.']
    try:
        refAF = alleles.count('0') / len(alleles)
        if minor_alt:
            maf = 1 - refAF
        else:
            maf = refAF
        maf = round(maf, 6)
    except ZeroDivisionError:
        maf = 'NA'
    #Count Genotypes
    genotypes = Counter(groupList)
    homRef = genotypes.get('0/0', 0)
    miss = genotypes.get('./.', genotypes.get('.',0))
    het = 0
    homVar = 0
    for genos, count in genotypes.items():
        if genos != '0/0' and genos != './.' and genos != '.':
            geno = genos.split('/')
            if geno[0] == geno[1]:
                homVar += count
            else:
                het += count
    
    return(','.join(map(str,[homRef, het, homVar, miss, maf])))

def createHeader(sortedGroups):
    first = "chrm,pos,filter,ref,alt,major,minor,minor_alt" 
    grp_suffix = ['count.homref', 'count.het', 'count.homvar', 'count.nocall', 'maf']
    
    middle = []
    for i in sortedGroups:
        for k in grp_suffix:
            middle.append(f"{i}.{k}")
    
    middle = ','.join(middle)
  
    end = (
    "Allele,Consequence,IMPACT,Gene,Feature_type,Feature,BIOTYPE,EXON,INTRON,cDNA_position,"
    "CDS_position,Protein_position,Amino_acids,Codons,STRAND,VARIANT_CLASS,protein_id"
    )
  
    header = f"{first},{middle},{end}"
    return(header)

#parse dog list  
dogs = {}
groups = set()
groups.add("all")
for i in sys.argv[2:]:
    group = os.path.basename(i).split('.')[0]
    groups.add(group)
    with open(i) as f: 
        for sample in f:
            sample = sample.strip()
            dogs[sample] = group

sortedGroups = list(groups)
sortedGroups.sort()

print(createHeader(sortedGroups))

#process vcf

vcfR = vcf.Reader(filename = sys.argv[1])
for line in vcfR:
    chrm, pos, ref, alt = line.CHROM, line.POS, line.REF, '/'.join(map(str, line.ALT))
    if line.FILTER:
        filt = line.FILTER[0]
    else:
        filt = 'PASS'
    AF = line.INFO.get('AF')
    AF = sum(AF)
    minor_alt = True
    if AF < 0.5:
        major, minor = ref, alt
    else:
        major, minor = alt, ref
        minor_alt = False 
    first = ','.join(map(str, [chrm, pos, filt, ref, alt, major, minor, minor_alt]))
    group_lists = {group: [] for group in groups}
    for i in line.samples:
        sample, gt = i.sample, i['GT'].replace('|', '/')
        group_lists['all'].append(gt)
        group = os.path.basename(dogs.get(sample, "")) if dogs.get(sample) is not None else None
        if group in group_lists:
            group_lists[group].append(gt) 
    groupAfs = {g: calcAF(gts, minor_alt) for g, gts in group_lists.items()}
    afs_str = ','.join([groupAfs[i] for i in  sortedGroups])
    csq = reduceCSQ(line.INFO.get('CSQ', []))
    for consequence in csq:
        print(f"{first},{afs_str},{consequence}")
