#! /usr/bin/env python

import sys
import pandas as pd

dat = pd.read_csv(sys.argv[1], sep = '\s+')

dat["het_rate"] = (dat["N(NM)"] - dat["O(HOM)"]) / dat["N(NM)"]

het_mean = dat["het_rate"].mean()
het_std = dat["het_rate"].std()

het_fail = dat[(dat["het_rate"] < het_mean - 3 * het_std) | (dat["het_rate"] > het_mean + 3 * het_std)]

het_fail["het_dst"] = (het_fail["het_rate"] - het_mean) / het_std

het_fail.to_csv("fail-het-qc.txt", sep="\t", index=False)
