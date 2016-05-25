#!/usr/bin/env python

import sys
import pandas as pd 

fn = sys.argv[1]
df = pd.read_table(fn)
out = open("../data/kmer/kebabs_kmer_list/"+fn.split("/")[-1].split(".")[0]+".list", "w+")

for kmer in df.index.values:
  out.write(kmer+"\n")

print "Writing to "+ "../data/kmer/kebabs_kmer_list/"+fn.split("/")[-1].split(".")[0]+".list"