#! /usr/bin/env python

import sys
import pandas as pd

df = pd.read_csv(sys.argv[1])

df.to_excel(sys.argv[2], index = False, engine='xlsxwriter')

