#! /usr/bin/env python

from cyvcf2 import VCF
import sys
import os
import argparse
from collections import defaultdict
from scipy.stats import fisher_exact


def fishers_test(geno1, geno2):
  grp1_ref, grp1_het, grp1_hom = geno1.count(0), geno1.count(1), geno1.count(2)
  grp2_ref, grp2_het, grp2_hom = geno2.count(0), geno2.count(1), geno2.count(2)

  f_rec_table = [[grp1_hom, grp1_het + grp1_ref],
                 [grp2_hom, grp2_het + grp2_ref]]

  f_dom_table = [[grp1_hom + grp1_het, grp1_ref],
                 [grp2_hom + grp2_het, grp2_ref]]
  
  if sum(geno1 + geno2) > 0:
    f_odds_rec, p_rec = fisher_exact(f_rec_table)
    f_odds_dom, p_dom = fisher_exact(f_dom_table)
  else:
    p_dom, p_rec = 1.0, 1.0

  return round(p_rec, 10), round(p_dom, 10)

def reduce_CSQ(csqs, csq_indices):
  if isinstance(csqs, str):
    csqs = csqs.split(',')

  updated = []
  for conseq in csqs:
    record = conseq.split('|')
    new_record = []
    for i in record:
      if '&' in i:
        new_record.append(i.split('&')[0])
      else:
        new_record.append(i)    
    result = ','.join([new_record[i] for i in csq_indices])
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

def fix_genos(genos, minor_alt, alt_cmd):

  new_geno = []
  if alt_cmd == True:
    for i in genos:
      first, second = (i[0], i[0]) if len(i) < 3 else (i[0], i[1])
      new_geno.append(-1 if -1 in (first, second) else
                        1 if first != second else
                        (0 if first == 0 else 2))

  else:
    for i in genos:
      first, second = (i[0], i[0]) if len(i) < 3 else (i[0], i[1])
      new_geno.append(-1 if -1 in (first, second) else 
                          1 if first != second else 
                          (0 if minor_alt else 2) if first == 0 else 
                          (2 if minor_alt else 0))
            
  return new_geno

def calc_af(genos, minor_alt, alt_cmd):
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
    if alt_cmd == True: 
      maf = 1 - refAF
    else:
      if minor_alt:
        maf = 1 - refAF
      else:
        maf = refAF
      maf = round(maf, 4)
  except ZeroDivisionError:
    maf = 'NA'

  return ','.join(map(str, [hom_ref, het, hom_var, missing, maf]))
 
def create_header(groups_sorted, csq_fields, target_sams, sv_file, alt_cmd, fishers):
  if alt_cmd:
    af = "alt_af"
  else:
    af = "minor_af"

  if sv_file:
    first = "chrm,pos,filter,ref,alt,major,minor,minor_alt,end,svtype,svlen"
  else:
    first = "chrm,pos,filter,ref,alt,major,minor,minor_alt"
  
  grp_suffix = ['count.homref', 'count.het', 'count.homvar', 'count.nocall', af]

  middle = []
  for i in groups_sorted:
    for k in grp_suffix:
      middle.append(f"{i}.{k}")
  
  middle = ','.join(middle)
  end = ','.join(csq_fields)
  if fishers:
    middle = f"{middle},fishers_p_rec,fishers_p_dom"

  sams = ','.join(target_sams)
  header = f"{first},{middle},{end},{sams}"

  return header

def create_groups_from_lists(direct_path):
  samples = defaultdict()
  groups = set()
  groups.add("all")
  for filename in os.listdir(direct_path):
    group_name = filename.split('.')[0]
    groups.add(group_name)
    filepath = os.path.join(direct_path, filename)
    if os.path.isfile(filepath):
      with open(filepath) as f:
        for sample in f:
          sample = sample.strip()
          samples[sample] = group_name
  
  sorted_groups = list(groups)
  sorted_groups.sort()

  return sorted_groups, samples

def check_csq(print_csqs, vcf_file, des_to_get):
  if print_csqs:
    for header in vcf_file.header_iter():
      if header.info().get('ID') == 'CSQ':
        csq_des = header.info().get("Description")
        if csq_des:
          desc, csq_fields = csq_des.split(':')
          csq_fields = csq_fields.replace(' ', '').split('|')
          lacking = list(set(des_to_get) - set(csq_fields))
          lacking = '\n'.join(lacking)
          csq_string = '\n'.join(csq_fields)
          print(f"Found CSQ attributes:\n{csq_string}\n")
          if lacking:
            print(f"Count not find required fields in CSQ attributes:\n{lacking}\nEdit des_to_get in script to match CSQ fields")
        else:
          print("Could not find VEP CSQ field in header")

    sys.exit()
  else:
    pass

def get_args():
  parser = argparse.ArgumentParser(description="counts the genotypes and allele frequency for defined groups")
  
  parser.add_argument('-v', '--vcf', required = True, help = "Path to VCF file")
  parser.add_argument('-d', '--dir', required = True, help = "Path to directory containing list files")
  parser.add_argument('-r', '--region', default = None, help = "Optional: region to analyze. Format must be chr or chr:start-end")
  parser.add_argument('-c', '--consequences', action = 'store_true', help = "Find and print all VEP CSQ attributes")
  parser.add_argument('-f', '--fishers', nargs = 2, metavar = ('GROUP1', 'GROUP2'), help = "Perform Fisher's Exact Test between two groups. Group names must be identical to file names in -d <dir>")
  parser.add_argument('-s', '--sv', action = 'store_true', help = "Use for SV vcf files")
  parser.add_argument('-a', '--alt', action = 'store_true', help = "Genotypes and allele frequencies are oriented toward minor allele. Add this arguement to orient toward alternate allele")  
  return parser.parse_args()

def main():
  #extrat argguments
  args = get_args()
  
  vcf = args.vcf
  list_dir = args.dir
  region = args.region
  print_csqs = args.consequences
  sv_file = args.sv
  alt_cmd = args.alt 
  fisher_groups = args.fishers

  #read vcf 
  vcf_file = VCF(vcf)
  vcf_samples = vcf_file.samples
  
  #get these main fields from csq  
  if sv_file:
    #SV VCF CSQ descriptions
    des_to_get = ['Allele', 'Consequence', 'IMPACT', 'Gene', 'Feature_type', 'Feature', 'BIOTYPE', 'EXON', 'INTRON']
  else:
    #SNP VCF CSQ descriptions
    des_to_get = ['Allele', 'Consequence', 'IMPACT', 'Gene', 'Feature_type', 'Feature', 'BIOTYPE', 'EXON', 'INTRON',
                'cDNA_position','CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'VARIANT_CLASS', 'ENSP', 
                'PhyloP_score', 'SIFT_SCORE', 'SIFT_PREDICTION']
  
  #If -c option, prints the CSQ fiels found in the VCF and missing CSQ feilds
  check_csq(print_csqs, vcf_file, des_to_get) 

  #gets inex of vcf file 
  csq_indices = extract_csq_indices(vcf_file, des_to_get)
  
  #creae and add samples to groups based on the list they're in
  groups_sorted, samples = create_groups_from_lists(list_dir)

  if fisher_groups: 
    #check for missing names if using a Fishers test 
    f_grp1, f_grp2 = fisher_groups
    missing = [g for g in (f_grp1, f_grp2) if g not in groups_sorted]
    if missing:
      sys.exit(f"Error: group name(s) not found: {', '.join(missing)}\nAvailable groups: {', '.join(groups_sorted)}") 

  gt_sample_order = [i for i in samples.keys() if i in vcf_samples]

  #create and print the header
  print(create_header(groups_sorted, des_to_get, gt_sample_order, sv_file, alt_cmd, fisher_groups))

  
  #parse through VCF
  for record in vcf_file(region):
    chrm, pos, ref, alt = record.CHROM, record.POS, record.REF, '/'.join(map(str,record.ALT))

    if record.FILTER:
      filt = record.FILTER
    else:
      filt = 'PASS'

    #AF need to be present in vcf, use bcftools +fill-tags to populate
    af = record.INFO.get('AF')
    if isinstance(af, tuple):
      af = sum(af) 
    minor_alt = True
    if af < 0.5:
      major, minor = ref, alt
    else:
      major, minor = alt, ref
      minor_alt = False

    if sv_file:
      svtype = record.INFO.get('SVTYPE')
      end = record.INFO.get('END')
      svlen = record.INFO.get('SVLEN') 
      first = ','.join(map(str, [chrm, pos, filt, ref, alt, major, minor, minor_alt, end, svtype, svlen]))
    else:
      first = ','.join(map(str, [chrm, pos, filt, ref, alt, major, minor, minor_alt]))

    group_lists = {pheno: [] for pheno in groups_sorted}
    header_order = []
    var_geno = dict(zip(vcf_samples, record.genotypes))

    for sample, gt in var_geno.items():
      group_lists['all'].append(gt)
      group_lists[samples.get(sample, "")].append(gt) if samples.get(sample) is not None else None
      if sample in gt_sample_order:
        header_order.append(gt)

    main_genos = fix_genos(header_order, minor_alt, alt_cmd)

    no_miss = [i for i in main_genos if i != -1]
    total_minor = sum(no_miss)
    
    fixed = ','.join(map(str, main_genos))

    group_afs = {g: calc_af(gts, minor_alt, alt_cmd) for g, gts in group_lists.items()} 
    afs_str = ','.join([group_afs[i] for i in groups_sorted])

    if fisher_groups:
      grp1_genos = fix_genos(group_lists[f_grp1], minor_alt, alt_cmd)
      grp2_genos = fix_genos(group_lists[f_grp2], minor_alt, alt_cmd) 
      f_rec_p, f_dom_p = fishers_test(grp1_genos, grp2_genos)
      afs_str = f"{afs_str},{f_rec_p},{f_dom_p}"
    
    csq = reduce_CSQ(record.INFO.get('CSQ', []), csq_indices)
    for cons in csq:
      print(f"{first},{afs_str},{cons},{fixed}")

if __name__ == "__main__":
  main() 
