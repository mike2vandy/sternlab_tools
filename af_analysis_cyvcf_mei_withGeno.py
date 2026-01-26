#! /usr/bin/env python

from cyvcf2 import VCF
import sys
import os
import argparse

def reduceCSQ(csqs, csq_indices):
  if isinstance(csqs, str):
    csqs = csqs.split(',')

  updated = []
  for conseq in csqs:
    record = conseq.split('|')
    result = ','.join([record[i] for i in csq_indices])
    updated.append(result)

  return updated

def extract_csq_indices(vcf_file, fields):
  for header in vcf_file.header_iter():
    if header.info().get('ID') == 'CSQ':
      csq_des = header.info().get("Description")
      break
  
  csq_fields = csq_des.split(" ")[-1].replace('"', '').split('|')
  csq_indices = [csq_fields.index(i) for i in fields]

  return csq_indices 

def fixGenos(genos, minor_alt):
  new_geno = []
  for i in genos:
    first, second = (i[0], i[0]) if len(i) < 3 else (i[0], i[1])
    new_geno.append(-1 if -1 in (first, second) else
                        1 if first != second else
                        (0 if minor_alt else 2) if first == 0 else
                        (2 if minor_alt else 0))

  return new_geno

def calcAF(genos, minor_alt):
  alleles = []
  hom_ref = 0
  het = 0
  hom_var = 0
  missing = 0
  for i in genos:
    if len(i) == 3:
      first, second, phased = i
    else:
      first = second = i[0]
    if first != -1 and second != -1:
      alleles.append(first)
      alleles.append(second)
      if first == 0 and first == second:
        hom_ref += 1
      elif first != second:
        het += 1
      else:
        hom_var += 1
    else:
      missing += 1

  try:
    refAF = alleles.count(0) / len(alleles)
    if minor_alt:
      maf = 1 - refAF
    else:
      maf = refAF
    maf = round(maf, 4)
  except ZeroDivisionError:
    maf = 'NA'

  return ','.join(map(str, [hom_ref, het, hom_var, missing, maf]))
 
def createHeader(sortedGroups, csq_fields, target_sams):
  first = "chrm,pos,filter,alt,minor_alt,svtype,mei,svlen"
  grp_suffix = ['count.homref', 'count.het', 'count.homvar', 'count.nocall', 'maf']

  middle = []
  for i in sortedGroups:
    for k in grp_suffix:
      middle.append(f"{i}.{k}")
  
  middle = ','.join(middle)
  end = ','.join(csq_fields)
  sams = ','.join(target_sams)
  header = f"{first},{middle},{end},{sams},total_minor"

  return header

def createGroupsFromLists(direct_path):
  samples = {}
  groups = set()
  groups.add("all")
  for filename in os.listdir(direct_path):
    group = filename.split('.')[0]
    groups.add(group)
    filepath = os.path.join(direct_path, filename)
    if os.path.isfile(filepath):
      with open(filepath) as f:
        for sample in f:
          sample = sample.strip()
          samples[sample] = group
  
  sortedGroups = list(groups)
  sortedGroups.sort()

  return sortedGroups, samples

def get_args():
  parser = argparse.ArgumentParser(description="counts the genotypes and allele frequency for groups")
  
  parser.add_argument('-v', '--vcf', required = True, help = "Path to VCF file")
  parser.add_argument('-d', '--dir', required = True, help = "Path to directory containing list files")
  parser.add_argument('-r', '--region', default = None, help = "Optional: region to analyze. Format must be chr or chr:start:end")
  
  return parser.parse_args()

def main():
  #extrat argguments
  args = get_args()
  
  vcf = args.vcf
  list_dir = args.dir
  region = args.region
  
  #rad vcf 
  vcf_file = VCF(vcf)
  vcf_samples = vcf_file.samples
  
  #get these main fields from csq 
  des_to_get = ['Allele', 'Consequence', 'IMPACT', 'Gene', 'Feature_type', 'Feature', 'BIOTYPE', 'EXON', 'INTRON']
    
  #get indices for CSQ descriptions
  csq_indices = extract_csq_indices(vcf_file, des_to_get)
  
  #creae and add samples to groups based on the list they're in
  sortedGroups, samples = createGroupsFromLists(list_dir)
  arry_sams = []
  for i in vcf_samples:
    if i in samples:
      if samples[i] == 'tvd':
        arry_sams.append(i)
  #create and print the header
  print(createHeader(sortedGroups, des_to_get, arry_sams)) 
  
  #parse through VCF
  for record in vcf_file(region):
    chrm, pos, alt = record.CHROM, record.POS, record.ID
    if 'abs' in alt:
      alt = 'DEL'
    elif 'ins' in alt:
      alt = 'INS'
    if record.FILTER:
      filt = record.FILTER
    else:
      filt = 'PASS'
    af = record.INFO.get('AF')
    svtype = record.INFO.get('SVTYPE')
    mei = record.INFO.get('MEI')
    svlen = record.INFO.get('SVLEN')
    if isinstance(af, tuple):
      af = sum(af)  
    minor_alt = True
    if af > 0.5:
      minor_alt = False
    first = ','.join(map(str, [chrm, pos, filt, alt, minor_alt, svtype, mei, svlen]))
    group_lists = {pheno: [] for pheno in sortedGroups}
    var_geno = dict(zip(vcf_samples, record.genotypes))
    for sample, gt in var_geno.items():
      group_lists['all'].append(gt)
      group_lists[samples.get(sample, "")].append(gt) if samples.get(sample) is not None else None
    shi_genos = group_lists['tvd']
    fixed = fixGenos(shi_genos, minor_alt)
    no_miss = [i for i in fixed if i != -1]
    total_minor = sum(no_miss)
    fixed = ','.join(map(str, fixed))
    groupAFs = {g: calcAF(gts, minor_alt) for g, gts in group_lists.items()}
    afs_str = ','.join([groupAFs[i] for i in sortedGroups])
    csq = reduceCSQ(record.INFO.get('CSQ', []), csq_indices)
    for cons in csq:
      print(f"{first},{afs_str},{cons},{fixed},{total_minor}")
  
if __name__ == "__main__":
  main() 
