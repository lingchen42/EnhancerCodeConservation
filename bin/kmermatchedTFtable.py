#!/usr/bin/env python
#function: get the top (FRAC) scored k-mers (FRAC*1024) and random (FRAC*1024) k-mers from the model, find the TFs they matched and put into a table.

import pandas as pd
import sys 
import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np
import os

def assign(row):
  #assign category
  if row["index"] <= cutoff_neg: 
    return "top_neg"
  elif row["index"] > cutoff_neg and row["score"] < 0:
    return "other_neg"
  elif row["index"] >= cutoff_pos:
    return "top_pos"
  elif row["index"] < cutoff_pos and row["score"] >= 0:
    return "other_pos"

def RC(seq):
  rc_dict = {}
  rc_dict["A"] = "T"
  rc_dict["C"] = "G"
  rc_dict["G"] = "C"
  rc_dict["T"] = "A"
  rc = []
  for i in seq:
    rc.append(rc_dict[i])
  return "".join(rc[::-1])

def tfs(category, df, rc=None):

  grped = df.groupby("category")
  grp = grped.get_group(category)
  
  tfs = []
  
  for value in grp["matched_TFs"].values:
    tfs = tfs + value

  temp = pd.DataFrame()

  if rc:
    temp[category+"rc"] = list(set(tfs))
    temp = pd.concat([DF_TF,temp],axis=1)
  else:
    temp[category] = list(set(tfs))
    temp = pd.concat([DF_TF,temp], axis=1)

  return temp


def rantfs(df,rc=None):
  temp = pd.DataFrame()

  for n in range(100):
    np.random.seed(n)
    rand_index = np.random.choice(len(df),int(FRAC*len(df)))

    tfs = []
    for value in df.ix[rand_index]["matched_TFs"].values:
      tfs = tfs + value


    temp2 = pd.DataFrame()
    if rc:
      temp2["random_rc%s"%(n)] = list(set(tfs))
      temp =pd.concat([temp,temp2], axis=1)

    else:
      temp2["random%s"%(n)]= list(set(tfs)) 
      temp =pd.concat([temp,temp2],axis=1)

  temp = pd.concat([DF_TF,temp],axis=1)
  return temp



MATCH_FN = sys.argv[1]
KMERSCORE_FN = sys.argv[2]
EXPRESSION_FN = sys.argv[3]
OUT_DIR = sys.argv[4]
if not os.path.exists(OUT_DIR):
    os.makedirs(OUT_DIR)
FRAC = 0.05


################################################################################
#kmer order
kmerorder_df = pd.read_table(KMERSCORE_FN)
kmers = kmerorder_df.index
kmer_score_dict = dict(zip(kmers, kmerorder_df["V1"]))

################################################################################
#liver expression dict
ex_df = pd.read_table(EXPRESSION_FN)
liver_dict =dict(zip(ex_df["gene_name"],ex_df["Liver"]))

################################################################################
#add jascore name dict 2014
TFname_dict={}
jascore_name = open("../data/motifdatabase/JASPAR_CORE_2014_vertebrates_list.txt", "r")
for line in jascore_name.readlines():
  eles = line[:-1].split(" ")
  TFname_dict[eles[0]] = eles[1]

################################################################################
# count the matches of highweightkmers, compared that to the rest. (highweigh, positive 20% vs postive 80% in ~ 512 positives; same for the negatives)
match_df = pd.read_table(MATCH_FN)
grped = match_df.groupby("#Query ID")

df_stats=pd.DataFrame(columns=["kmer", "score","matched_TFs", "matched_TFs_counts","liver_expression", "liver_expressed_TF_counts"])
for kmer in kmers:
  try:
    grp = grped.get_group(kmer)
    matchedTFs = []
    for TF in grp["Target ID"].values:
      try:
        matchedTFs.append(TFname_dict[TF])
      except KeyError:
        matchedTFs.append(TF.upper())

    liver_expression = []
    for TF in matchedTFs:
      try:
        liver_expression.append(liver_dict[TF])
      except KeyError:
        liver_expression.append("NA")
    liver_c = [x for x in liver_expression if x not in [0, "NA"] ]

  except KeyError:
    matchedTFs = []
    liver_expression = []
    liver_c = []

  score = kmer_score_dict[kmer]

  row = [kmer, score, matchedTFs, len(matchedTFs), liver_expression, len(liver_c)]
  df_stats.loc[len(df_stats)+1] = row

cutoff_neg = int(len(df_stats)*FRAC)
cutoff_pos = len(df_stats) - cutoff_neg
df_stats = df_stats.sort("score")
df_stats.reset_index(inplace=True)
df_stats["category"] = df_stats.apply(assign,axis=1)

df_stats.to_csv(OUT_DIR+MATCH_FN.split("/")[-2].replace(".meme","_matchedTF.csv"))
print "Writing "+OUT_DIR+MATCH_FN.split("/")[-2].replace(".meme","_matchedTF.csv")

################################################################################
#consider RC
df_stats_rc=pd.DataFrame(columns=["kmer","score","matched_TFs","matched_TFs_counts","liver_expression","liver_expressed_TF_counts"])
kmer_rc = []
for kmer in kmers:
  if RC(kmer) in kmer_rc:
    continue
  else:
    kmer_rc.append(kmer)

kmer_match_dict = dict(zip(df_stats["kmer"], zip(df_stats["score"],df_stats["matched_TFs"],df_stats["matched_TFs_counts"], df_stats["liver_expression"],df_stats["liver_expressed_TF_counts"])))
for kmer in kmer_rc:

  score = kmer_match_dict[kmer][0]+kmer_match_dict[RC(kmer)][0]
  matchedTFs = kmer_match_dict[kmer][1]+kmer_match_dict[RC(kmer)][1]
  matchedTFs_c = kmer_match_dict[kmer][2]+kmer_match_dict[RC(kmer)][2]
  liver = kmer_match_dict[kmer][3]+kmer_match_dict[RC(kmer)][3]
  liver_c = kmer_match_dict[kmer][4]+kmer_match_dict[RC(kmer)][4]
  row = [kmer, score, matchedTFs, matchedTFs_c, liver, liver_c]
  df_stats_rc.loc[len(df_stats_rc)+1] = row


df_stats_rc = df_stats_rc.sort("score")
df_stats_rc.reset_index(inplace=True)
cutoff_neg = int(len(df_stats_rc)*FRAC)
cutoff_pos = len(df_stats_rc) - cutoff_neg
df_stats_rc["category"] = df_stats_rc.apply(assign,axis=1)


################################################################################
#build a TF table
#TOP_POS_MATCHED_TFS TOP_NEG_MATCHED_TFS RANDON_KMER_MATCHED_TFS_1 RANDON_KMER_MATCHED_TFS_2 RANDON_KMER_MATCHED_TFS_3 .... TOP_POS_MATCHED_TFS_RC TOP_NEG_MATCHED_TFS_RC RANDON_KMER_MATCHED_TFS_RC

DF_TF = pd.DataFrame()

categories = ["top_pos","top_neg"]
for category in categories:
  DF_TF=tfs(category,df_stats)
for category in categories:
  DF_TF=tfs(category,df_stats_rc,rc=1)

DF_TF=rantfs(df_stats)
DF_TF=rantfs(df_stats_rc,rc=1)

DF_TF.to_csv(OUT_DIR+KMERSCORE_FN.split("/")[-1].split(".")[0]+"_TFmatches.tab", sep="\t", index=False)
print "Writing to ", OUT_DIR+KMERSCORE_FN.split("/")[-1].split(".")[0]+"_TFmatches.tab"







