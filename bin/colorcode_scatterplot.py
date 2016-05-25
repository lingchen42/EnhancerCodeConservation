#!/usr/bin/env python
# Draw a scatter plot with the SVM score of every enhancer, color code them with different properties.
# usage ./colorcode_scatterplot.py -fi /dors/capra_lab/data/dna/hg19.fa -sf K5/Hsap_cv_scores_2015-08-07.tsv K5/Btau_applyto_Hsap_2015-08-27.tsv  -r ../randomGenomicRegion/Hsap/Hsap_Enhancers_gapexcluded.bed4 -c gc -con ../../villar15/highlyConserved/highlyConservedEnhancers_H3K27ac_placentals ../../villar15/lineageSpecific/Hsap_LineageSpEnhancers_H3K27ac_primates_epo ../../villar15/recentlyEvolved/Hsap_H3K27ac_humanspEnhancers 

import os
import datetime
#import pybedtools
import argparse
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
import pandas as pd
from decimal import Decimal


#parse args
arg_parser = argparse.ArgumentParser(description="Colorcoded Scatterplot of SVM enhancer scores")
arg_parser.add_argument("-o", "--out_dir",help="output directory default=../figuredraft/enhancerscorecor/")
arg_parser.add_argument("-fi", "--fasta_file", help="fasta file to get the sequence")
arg_parser.add_argument("-gc",help="GC content table")
arg_parser.add_argument("-ph",help="H3K27ac peak height file")
arg_parser.add_argument("-sf", "--svm_score_file", nargs=2, help="a pair of score files of same sepcies enhancers in different models")
arg_parser.add_argument("-con","--conservation", nargs=3, help="original files of highlyConserved enhancers, lineageSpecific and recentlyEvolved.")
arg_parser.add_argument("-r","--region_file", help="bed file of the regions want to be plotted")
arg_parser.add_argument("-c","--color_code",choices =["gc","length","conservation", "bin_length","ph"],help="choose from gc, length, conservation, and peak height, which property to be color coded. bin_length is only available for table inputs")
arg_parser.add_argument("-t","--from_table", nargs=2,help="if the tables with other properties exit, input the tables")
args = arg_parser.parse_args()

#set global values
OUT_DIR = "../results/figures/scatterplots/"
if args.out_dir:
    OUT_DIR = os.path.join(args.out_dir, str(datetime.datetime.now())[:10])+'/'
    if not os.path.exists(OUT_DIR):os.makedirs(OUT_DIR)
if args.fasta_file:
  FASTA = args.fasta_file
else:
  FASTA = None
FHS = args.svm_score_file
if args.gc:
  GC = args.gc
else:
  GC = None
if args.conservation:
  CON = args.conservation
else:
  CON = None
if args.ph:
  PH = args.ph
else:
  PH = None
if args.region_file:
  REGIONS = args.region_file
else:
  REGIONS = None
if args.from_table:
  TABLE = args.from_table
else:
  TABLE = None

CHOICE = args.color_code


############################################################################
def name(row):
  name = "chr"+row["Chrom"]+":"+str(row["Start"])+"-"+str(row["End"])
  return name

def gc_frac(seq):    
  num_bases = float(len(seq) - seq.count('N'))
  return (seq.count('C') + seq.count('G')) / num_bases if num_bases > 0 else None

def read_con(conservation_files):
  #input the highlyConserved, lineageSepcific and recentlyEvolved in sequence
  hc_regions=[]
  ls_regions=[]
  re_regions=[]
  for index, fn in enumerate(CON):
    fh = open(fn, "r")
    for line in fh.readlines():
      if line.startswith("Chrom"):continue
      t = line.split("\t")
      if index == 0: hc_regions.append(["chr"+t[0]+":"+t[1]+"-"+t[2]])
      elif index == 1: ls_regions.append("chr"+t[0]+":"+t[1]+"-"+t[2])
      else: re_regions.append("chr"+t[0]+":"+t[1]+"-"+t[2]) 
  return hc_regions, ls_regions, re_regions

def bin_length(row):
  length = row["length"]
  #get the length bin
  if length <= 2000:
    return 0
  elif length > 2000 and length <= 4000:
    return 1
  elif length > 4000 and length <= 6000:
    return 2
  elif length > 6000 and length <= 8000:
    return 3
  else:
    return 4
 
def colorcode_scatter(df_1, df_2):
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

  if CHOICE == "conservation":
    dict_1 = dict(zip(df_1["names"], zip(df_1["decision"],df_1["conservation"])))
    dict_2 = dict(zip(df_2["names"], df_2["decision"]))

  elif CHOICE == "gc":
    dict_1 = dict(zip(df_1["names"], zip(df_1["decision"],df_1["gc"])))
    dict_2 = dict(zip(df_2["names"], df_2["decision"]))

  elif CHOICE == "length":
    dict_1 = dict(zip(df_1["names"], zip(df_1["decision"],df_1["length"])))
    dict_2 = dict(zip(df_2["names"],df_2["decision"]))
  elif CHOICE == "ph":
    dict_1 = dict(zip(df_1["names"], zip(df_1["decision"],df_1["peakheight"])))
    dict_2 = dict(zip(df_2["names"],df_2["decision"]))
  else:
    dict_1 = dict(zip(df_1["names"], zip(df_1["decision"],df_1["bin_length"])))
    dict_2 = dict(zip(df_2["names"],df_2["decision"]))


  for enhancer in region_names:
    try:
      s2_decision.append(dict_2[enhancer])
      s1_decision.append(dict_1[enhancer][0])
      s1_other.append(dict_1[enhancer][1])
    except KeyError:
      continue
      #s2_decision.append(np.nan)

  r_value, p_value = stats.spearmanr(s1_decision,s2_decision)

  #x, y labels are the names of species 
  n1 = FHS[0].split("/")[-1].split(".")[0]
  n2 = FHS[1].split("/")[-1].split(".")[0]

  s1 = n1.split("_")[0]
  s2 = n2.split("_")[0]

  fig = plt.figure()
  if s2=="Mmul":
    fig.set_size_inches(10, 8)
  else:
    fig.set_size_inches(8, 8)
  ax = fig.add_subplot(111)
  scatter = ax.scatter(s1_decision, s2_decision, c= s1_other, alpha=0.5)

  plt.xlabel(s1+" classifier")
  plt.ylabel(s2+" classifier")
  
  if n1.find("gc")==-1:
    plt.xlim(-6,6)
    plt.ylim(-6,6)

  plt.text(0.80,0.95, r"$\rho =%s$"%("{0:.2f}".format(r_value)), ha='center', va='center', transform=ax.transAxes, fontsize=25)

  if (s2=="Mmul")&(CHOICE=="gc"):
    cbar = plt.colorbar(scatter, label="GC")
  if (s2=="Mmul")&(CHOICE=="length"):
    #colorbar for binned length
    bounds = [0,1,2,3,4,5]
    labels = [0,2000,4000,6000,8000,""]
    cbar = plt.colorbar(scatter, boundaries=bounds, ticks=bounds, label="Length")
    cbar.set_ticklabels(labels)

  matplotlib.rcParams.update({'font.size': 22})

  if n1.find("gc")!=-1:
    plt.savefig(OUT_DIR+"%s_%s_%s_GCaware.pdf"%(s1,s2,CHOICE))
    print "Write to "+OUT_DIR+"%s_%s_%s_GCaware.pdf"%(s1,s2,CHOICE)
  else:
    plt.savefig(OUT_DIR+"%s_%s_%s.pdf"%(s1,s2,CHOICE))
    print "Write to "+OUT_DIR+"%s_%s_%s.pdf"%(s1,s2,CHOICE)

############################################################################

if not TABLE:
  ##First, make a dataframe which store all the properties
  ##Get GC and length for each region

  if FASTA:
    # fasta file parser
    fasta = pybedtools.example_filename(FASTA)
  
  if GC:
    df_gc = pd.read_table(GC)
    gc_dict = dict(zip(df_gc["names"],df_gc["gc"]))

  #conservation files
  #make 3 list have different conservation status: highly conserved, lineage specific and recently evovled
  if CON:
    hc_regions, ls_regions, re_regions = read_con(CON)

  #H3K27ac peak height info
  if PH:
    df_ph = pd.read_table(PH)
    df_ph["names"] = df_ph.apply(name,axis=1)
    ph_dict = dict(zip(df_ph["names"],df_ph["MeanFoldChange"]))
  
  #read in the SVM score files
  for index1, fh in enumerate(FHS):
  
    df = pd.read_table(fh)
    #Compute GC, length and conservation through iteration of rows
    lengths = []

    if FASTA or GC:
      gcs = []
    if CON:
      conservation = []
    if PH:
      phs = []
    
    for index2, row in df.iterrows():
      # get the coordinate
      region = row[0]
      t = region.split(":")
      chrom = t[0]
      start = int(t[1].split("-")[0])
      end = int(t[1].split("-")[1])

      # get the length
      length = end-start
      lengths.append(length)

      #get sequence through pybedtools
      if FASTA:
        if not start == end:
          region_bed = pybedtools.BedTool(" ".join([chrom, str(start), str(end)]), from_string=True)
          t_seq = region_bed.sequence(fi=fasta)
          region_seq = open(t_seq.seqfn).read().split("\n")[1]
          gcs.append(gc_frac(region_seq.upper()))
        else:
          gcs.append(np.nan)

      #get gc from table
      if GC:
        try:
          gcs.append(gc_dict[row[0]])	     
        except KeyError:
          gcs.append(np.nan)

      # get the conservation status
      if CON:
        #highlyConserved as 3
        if region in hc_regions: 
          conservation.append("3")
        #lineageSpecific as 2
        elif region in ls_regions: 
          conservation.append("2")
        #recentlyEvolved as 1
        elif region in re_regions: 
          conservation.append("1")
        #others as 0
        else:
          conservation.append("0")

      #get peak height
      if PH:
        try:
          phs.append(ph_dict[row[0]])
        except KeyError:
          phs.append(np.nan)


    df["length"] = lengths
    if FASTA or GC:
      df["gc"] = gcs
    if CON:
      df["conservation"] = conservation
    if PH:
      df["peakheight"] = phs

    if index1 == 0 : 
      df_1 = df
      #df_1.to_csv(OUT_DIR+fh.split("/")[-1].split(".")[0]+"_withotherfeature.tsv",sep="\t")
      df_1.to_csv("../results/SVM_output/withotherfeature/"+fh.split("/")[-1].split(".")[0]+"_withotherfeature.tsv",sep="\t")
      print "Writing to ../results/SVM_output/withotherfeature/"+fh.split("/")[-1].split(".")[0]+"_withotherfeature.tsv"
    elif index1 == 1:
      df_2 = df
      #df_2.to_csv(OUT_DIR+fh.split("/")[-1].split(".")[0]+"_withotherfeature.tsv",sep="\t")
      df_2.to_csv("../results/SVM_output/withotherfeature/"+fh.split("/")[-1].split(".")[0]+"_withotherfeature.tsv",sep="\t")
      print "Writing to ../results/SVM_output/withotherfeature/"+fh.split("/")[-1].split(".")[0]+"_withotherfeature.tsv"


  colorcode_scatter(df_1,df_2)

else:
  df_1 = pd.read_table(TABLE[0])
  df_2 = pd.read_table(TABLE[1])
  FHS = TABLE
  df_1["bin_length"] = df_1.apply(bin_length,axis=1)
  df_2["bin_length"] = df_2.apply(bin_length,axis=1)
  colorcode_scatter(df_1,df_2)





  

   






