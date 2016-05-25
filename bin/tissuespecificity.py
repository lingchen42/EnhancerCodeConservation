#!/usr/bin/env python
# input a list of TF, output the specificity score from http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2836267/#SD6
# widespread tissue expression (TSPS < 1) and a second smaller population at higher tissue specificity (TSPS >= 1)

import sys
import pandas as pd
import numpy as np
import os

tf_fh = open(sys.argv[1],"r")
tfs=[]
for line in tf_fh.readlines():
  tfs.append(line[:-1])

df_tf=pd.read_excel('../data/TFexpression/NIHMS177825-supplement-02.xls', header=1)

#parse Tf_with_facets.csv
genes=df_tf['Symbol (Human)'].values
spec=df_tf['Tissue specificity score'].values
#two lists to dictionary
#gene: facets
spec_dict=dict(zip(genes, spec))

ave = []

for tf in tfs:
  try:
    s = spec_dict[tf]
    if not pd.isnull(s):
      ave.append(s)
      print tf,"\t",s
  except KeyError:
    continue


## If we consider compound TFs as well. 
# for tf in tfs:
#   if tf.find("::")==1:
#     try:
#       s = spec_dict[tf]
#       if not pd.isnull(s):
#         ave.append(s)
#         print tf,"\t",s
#     except KeyError:
#       continue
#   else:
#     #for compound TFs, get their mean TSPS
#     eles= tf.split("::")
#     s=[]
#     for i in eles:
#       try:
#         s_i = spec_dict[tf]
#       except KeyError:
#         s_i.append(np.nan)
    
#     if not pd.isnull(np.mean(s)):
#       ave.append(s)
        

#print "mean specificity of %s TFs: "%(len(ave))
#print np.mean(ave),"+/-",np.std(ave)
