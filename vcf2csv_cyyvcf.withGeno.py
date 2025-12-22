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

def createHeader(csq_fields, target_sams):
  first = "chrm,pos,filter,ref,alt"
  end = ','.join(csq_fields)
  sams = ','.join(target_sams)
  header = f"{first},{end},{sams}"

  return header

def get_args():
  parser = argparse.ArgumentParser(description="counts the genotypes and allele frequency for groups")
  
  parser.add_argument('-v', '--vcf', required = True, help = "Path to VCF file")
  
  return parser.parse_args()

def main():
  #extrat argguments
  args = get_args()
  
  vcf = args.vcf 
  #rad vcf 
  vcf_file = VCF(vcf)
  vcf_samples = vcf_file.samples
  
  #get these main fields from csq 
  des_to_get = ['Consequence', 'IMPACT', 'Gene', 'PhyloP_score', 'SIFT_SIFT_SCORE', 'SIFT_SIFT_PREDICTION']
    
  #get indices for CSQ descriptions
  csq_indices = extract_csq_indices(vcf_file, des_to_get)
  
  #create and print the header
  print(createHeader(des_to_get, vcf_samples))
  
  #parse through VCF
  for record in vcf_file:
    chrm, pos, ref, alt = record.CHROM, record.POS, record.REF, '/'.join(map(str,record.ALT))
    if record.FILTER:
      filt = record.FILTER
    else:
      filt = 'PASS'
    first = ','.join(map(str, [chrm, pos, filt, ref, alt]))
    gts = ','.join(map(str,record.gt_types))
    gts = gts.replace('2','-1')
    gts = gts.replace('3','2')
    csq = reduceCSQ(record.INFO.get('CSQ', []), csq_indices)
    for cons in csq:
      print(f"{first},{cons},{gts}")

if __name__ == "__main__":
  main() 
