#!/usr/bin/env python

import os
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def parsespecies(fh):
  s = fh.split("/")[-1].split("_")[0]
  return s

OUT_DIR = "../results/figures/"
if not os.path.exists(OUT_DIR):
    os.makedirs(OUT_DIR)
outfile_name = OUT_DIR+"gcdist.pdf"

color_dict = {"Hsap":"#006600", "Mmul":"#5A9BD4","Mmus":"#F15A60","Btau":"#FFDC5B","Cfam":"#FF7500","Mdom":"#C09953"}

fhs = sys.argv[1:]
dfs = []
species = []
for fh in fhs:
  species.append(parsespecies(fh))
  dfs.append(pd.read_table(fh))

sns.set(font_scale=2.5)
sns.set_style("white")

fig = plt.figure()
fig.set_size_inches(8, 8)
ax1 = fig.add_subplot(111)

plt.xlim(0,1)

for i, df in enumerate(dfs):
  if species[i] == "Hsap":
    sns.distplot(df["gc"], label = species[i], color=color_dict[species[i]])
  sns.distplot(df["gc"], hist=False, label = species[i], color=color_dict[species[i]])


#plt.yticks(fig.get_yticks(), fig.get_yticks()/100)
plt.xlabel("GC")
plt.ylabel("Density")
plt.legend()

#matplotlib.rcParams.update({'font.size': 25})

plt.savefig(outfile_name)
print "Drawing ", outfile_name
