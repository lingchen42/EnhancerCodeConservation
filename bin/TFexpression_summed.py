#!/usr/bin/env python
# function <binary expression file> <selected_TF_list_fn> <all TF list>
import pandas as pd
import sys
from scipy import stats
import numpy as np

#binary expression profile
df = pd.read_table(sys.argv[1],sep=",")

#all TF with pwm
TF_with_pwm_fh = open(sys.argv[3], "r")
allTFwithpwm = []
for line in TF_with_pwm_fh.readlines():
  allTFwithpwm.append(line[:-1])

#customed TF list
TFlist_fn = open(sys.argv[2])
TFs = []
for line in TFlist_fn.readlines():
  TFs.append(line[:-1])

notsharedTFs = set(allTFwithpwm).difference(set(TFs))

def sumtissue(l):
  dft = df[df["gene_name"].isin(l)]
  x = dft.sum(axis=1,numeric_only=True)
  dft= dft.append(dft.sum(numeric_only=True), ignore_index=True)
  dft_T= dft.drop("gene_name",1).T

  print "Number of TFs: ", len(dft)-1
  print dft_T[dft_T.columns[-1]][:-1]
  return x

def sumtissue2(l):
  #consider compound TF motifs. If all tfs in the compound motifs have x-tissue expression, then the compound will be considered as x-tissue expressed.

  l1=[] #single TF motif
  l2=[] #compound TF motif

  for i in l:
    if i.find("::")==-1:
      l1.append(i)
    else:
      l2.append(i)

  dft = df[df["gene_name"].isin(l1)]

  for i in l2:
    eles = i.split("::")
    df_eles = df[df["gene_name"].isin(eles)]
    if len(df_eles)==len(eles):
      df_eles=df_eles.append(df_eles.sum(numeric_only=True), ignore_index=True).set_index("gene_name")[-1:]
      df_eles[df_eles!=len(eles)]=0
      df_eles[df_eles==len(eles)]=1
      df_eles = df_eles.reset_index()
      df_eles["gene_name"]=i
      dft=dft.append(df_eles,ignore_index=True)

  x = dft.sum(axis=1,numeric_only=True)
  dft= dft.append(dft.sum(numeric_only=True), ignore_index=True)
  dft_T= dft.drop("gene_name",1).T

  print "Number of TFs: ", len(dft)-1
  #print dft
  print dft_T[dft_T.columns[-1]][:-1]
  return x, dft


# print sys.argv[3].replace(".list", "")
# print "TFs not in input:"
# a = sumtissue(list(notsharedTFs))
# print "Input TFs:"
# b = sumtissue(TFs)

print sys.argv[3]
print "TFs not in input:"
a, df1 = sumtissue2(list(notsharedTFs))
print "Input TFs:"
b, df2 = sumtissue2(TFs)
out = sys.argv[2].replace(".list","_expressionprofile.csv")
df2.to_csv(out)
print "Writing Expression profile of input TFs into ", out


statistic, pvalue = stats.mannwhitneyu(a,b)
print "\n Mannwhitneyu test of the number of active tissues in two sets of TFs"
print "TFs not in in put active tissues mean: ", np.mean(a), "input TFs active tissues mean:", np.mean(b)
print "mannwhitneyu statistic: ",statistic, "p-value: ", pvalue