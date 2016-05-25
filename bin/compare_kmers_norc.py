#!/usr/bin/env python

import sys
import pandas as pd
import numpy as np

PERCENT = float(sys.argv[1])
fns = sys.argv[2:]

top_dict = {}
neg_tops = []
pos_tops = []

for fn in  fns:
  species = fn[:4]
  outname = fn.split("/")[-1].split(".")[0]
  df = pd.read_table(fn)
  num = int(PERCENT*len(df))

  if outname.find("Limb")!=-1:
    OUTDIR = "../results/figures/venn/topkmerlist_limb/"
  else:
    OUTDIR = "../results/figures/venn/topkmerlist/"

  neg_top = df.index[:num]
  neg_tops.append(set(neg_top))
  out = open(OUTDIR+outname+"_top%s_neg.list"%(sys.argv[1]),"w+")
  for i in neg_top:
    out.write(i+"\n")

  pos_top = df.index [-num:]
  pos_tops.append(set(pos_top))
  out = open(OUTDIR+outname+"_top%s_pos.list"%(sys.argv[1]),"w+")
  for i in pos_top:
    out.write(i+"\n")

  top_dict[species] = [pos_top, neg_top]


#shared among all species, rc counted
pos_all = set.intersection(*pos_tops)
neg_all = set.intersection(*neg_tops)
print "pos shared in all, no rc: ", pos_all
print "neg shared in all, no rc:", neg_all

species_list = top_dict.keys()
for i in species_list:
  for j in species_list:
    if species_list.index(i) < species_list.index(j): 
      i_pos = set(top_dict[i][0])
      j_pos = set(top_dict[j][0])
      pos_pair = i_pos&j_pos

      i_neg = set(top_dict[i][1])
      j_neg = set(top_dict[j][1])
      neg_pair = i_neg&j_neg

      middle = [len(pos_pair), len(neg_pair)]
      names = ["Top Pos","Top Neg"]
      df = pd.DataFrame()
      df["names"] = names
      #df["total"] = [num, num]
      df["shared in pair"] = [len(pos_pair), len(neg_pair)]
      df["shared in all"] = [len(pos_all), len(pos_all)]

      print i, " ", j
      print df, "\n"



