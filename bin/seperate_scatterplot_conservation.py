#!/usr/bin/env python

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy import stats
from scipy.stats import ttest_ind
import pandas as pd

FH1 = sys.argv[1]
FH2 = sys.argv[2]
REGIONS = sys.argv[3]
OUT_DIR = "../results/figures/scatterplots/"


def colorcode_scatter(df_1, df_2):
  n1 = FH1.split("/")[-1].split(".")[0]
  n2 = FH2.split("/")[-1].split(".")[0]

  s1 = n1.split("_")[0]
  s2 = n2.split("_")[0]

  #read in the regions need to be plotted
  if REGIONS:
    region_fh = open(REGIONS,"r")
    region_names = []
    for line in region_fh.readlines():
      eles=line.split("\t")
      region_names.append(eles[0]+":"+eles[1]+"-"+eles[2]) 
  else:
    region_names = df_1["names"]

  #pick one property to be color coded
  s1_decision = []
  s2_decision = []
  s1_other = []

  dict_1 = dict(zip(df_1["names"], zip(df_1["decision"],df_1["conservation"])))
  dict_2 = dict(zip(df_2["names"], df_2["decision"]))

  for enhancer in region_names:
    s1_decision.append(dict_1[enhancer][0])
    s1_other.append(dict_1[enhancer][1])
    try:
      s2_decision.append(dict_2[enhancer])
    except KeyError:
      s2_decision.append(np.nan)

  df_con = pd.DataFrame()
  df_con["enhancer_name"] = region_names
  df_con[n1+"_decision"] = s1_decision
  df_con[n2+"_decision"] = s2_decision
  df_con["conservation"] = s1_other

  df_con=df_con.dropna(how="any")
  grped = df_con.groupby("conservation")
  fig, axs = plt.subplots(3, 1, sharex=True, sharey=True)
  fig.set_size_inches(8, 24)
  if n1.find("gc")==-1:
    plt.xlim(-6,6)
    plt.ylim(-6,6)
  plt.xlabel(s1+" Classifier")
  axs[1].set_ylabel(s2+" Classifier")
  #plt.title("Enhancer Score Correlation")
  i=0
  color =["blue","yellow","red"]
  for key, grp in grped:
    slope, intercept, r_value, p_value, std_err = stats.linregress(grp[n1+"_decision"],grp[n2+"_decision"])
    
    axs[i].scatter(grp[n1+"_decision"], grp[n2+"_decision"], c= color[i])
    #axs[i].text(1,1,"r-squared %s"%(r_value), ha='right', va='top')
    i=i+1
  
  matplotlib.rcParams.update({'font.size': 25})

  if n1.find("gc")!=-1:
    plt.savefig(OUT_DIR+"%s_%s_%s_GCaware.pdf"%(n1,n2,"sep_conservation"))
    print "save to "+ OUT_DIR+"%s_%s_%s_GCaware.pdf"%(n1,n2,"sep_conservation")
  else:
    plt.savefig(OUT_DIR+"%s_%s_%s.pdf"%(n1,n2,"sep_conservation"))
    print "save to "+ OUT_DIR+"%s_%s_%s.pdf"%(n1,n2,"sep_conservation")

df_1 = pd.read_table(FH1)
df_2 = pd.read_table(FH2)
colorcode_scatter(df_1,df_2)