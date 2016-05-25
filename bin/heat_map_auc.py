#!/usr/bin/env python

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import sys

OUTDIR= "../results/figures/heatmaps/"
fh = sys.argv[1]
eles = fh.split("/")[-1].split(".")

if eles[-1] == "csv":
  df = pd.read_table(fh,sep=",",index_col=0)
else:
  df = pd.read_table(fh,sep="\t",index_col=0)

df = df.dropna(axis=0, how="all")
df = df.dropna(axis=1, how="all")

fig = plt.figure()
fig.set_size_inches(10, 8)
ax1 = fig.add_subplot(111)

if fh.find("relative") != -1:
  sns.heatmap(df, square=True, annot=True, annot_kws={"size": 20})
elif fh.find("decrease")!= -1:
  sns.heatmap(df, square=True, annot=True, cmap = "winter", annot_kws={"size": 20})
else:
  sns.heatmap(df, square=True, annot=True, cmap = "summer_r", annot_kws={"size": 20})

ax1.tick_params(labelsize=25)
cax = plt.gcf().axes[-1]
cax.tick_params(labelsize=25)

#print "Table mean value: ", np.mean(df.mean(numeric_only=True))

plt.savefig(OUTDIR+eles[0]+"_heatmap.pdf")
print "writing to ", OUTDIR+eles[0]+"_heatmap.pdf"

