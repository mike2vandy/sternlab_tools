#! /usr/bin/env python

import sys, gzip

header = '\n'.join(["##fileformat=VCFv4.2",
"##ALT=<ID=NON_REF,Description=\"Represents any possible alternative allele\">",
"##INFO=<ID=TRANSCRIPT,Description=\"Transcript ID from SIFT\">",
"##INFO=<ID=SIFT_SCORE,Description=\"SIFT score from SIFT\">",
"##INFO=<ID=SIFT_PREDICTION,Description=\"SIFT prediction from SIFT\">",
"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"])

print(header)
with gzip.open(sys.argv[1], 'rt') as f:
    for line in f:
        if line.startswith('#'): 
            continue
        else:
            chrm, pos, ref, alt, trans, score, pred = line.strip().split('\t')
            pred = pred.replace('low confidence', 'LC').replace(' ','')
            print(f"{chrm}\t{pos}\t.\t{ref}\t{alt}\t.\t.\tTRANSCRIPT={trans};SIFT_SCORE={score};SIFT_PREDICTION={pred}")
