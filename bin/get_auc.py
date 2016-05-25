#!/usr/bin/env python

import os
import subprocess
import pandas as pd 

def dict2tab(auc_dict,gc="",limb=False):
  
  tabs = []
  for i in range(4):
    tabs.append(pd.DataFrame(index=species, columns=species))

  for k in auc_dict:
    #df_roc
    tabs[0][k[0]][k[1]] = auc_dict[k][0] 
    #df_roc_relative
    tabs[1][k[0]][k[1]] = auc_dict[k][0] / auc_dict[(k[1],k[1])][0]
    # df_pr
    tabs[2][k[0]][k[1]] = auc_dict[k][1] 
    #df_pr_relative
    tabs[3][k[0]][k[1]] = auc_dict[k][1] / auc_dict[(k[1],k[1])][1] 

  names = ["roc","relative_roc","pr", "relative_pr"]
  for index, t in enumerate(tabs):

    if not limb:
      out_name = OUTDIR+names[index]+"%s.tab"%(gc)
    else:
      out_name = OUTDIR+"Limb_"+names[index]+"%s.tab"%(gc)
    
    t.to_csv(out_name, sep="\t")
    print "Writing to ", out_name


def getaucs(cmd,gc="",limb=False):

  output = subprocess.check_output(cmd, shell=True)

  t = output.split("\n\n")

  #(train, test): (roc_auc, pr_auc)
  auc_dict = {}

  for i in t[:-1]:
    i = [x for x in i.split("\n") if x]
    j= i[0].split("/")[-1].split("_")
    if limb:
      if i[0].find("cv") == -1:
        s1 = j[1]
        s2 = j[3]
      else:
        s1= j[1]
        s2= j[1]
    else:
      if i[0].find("cv") == -1:
        s1 =j[0]
        s2 =j[2]
      else:
        s1= j[0]
        s2= j[0]
    
    roc = float(i[1].split("ROC AUC:")[-1])
    pr = float(i[2].split("PR AUC:")[-1])
    auc_dict[(s1,s2)]=(roc, pr)

  dict2tab(auc_dict,gc,limb)

OUTDIR ="../results/figures/heatmaps/"
species = ["Hsap","Mmul","Mmus","Btau","Cfam","Mdom"]

files = [f for f in os.listdir('../results/SVM_output/') if f.find(".tsv") != -1 and f.find("Limb")==-1 and f.find("nommus")==-1 and f.find("without") == -1 and f.find("highweightkmers") == -1]
files_nogc = ["../results/SVM_output/"+f for f in files if f.find("gc") == -1 ]
files_gc = ["../results/SVM_output/"+f for f in files if f.find("gc") != -1 ]

cmd_nogc = " ".join(["./preds2roc_pr_curves_customed.py","t"] + files_nogc)
cmd_gc = " ".join(["./preds2roc_pr_curves_customed.py","t"] + files_gc)

getaucs(cmd_nogc)
getaucs(cmd_gc,gc="_gc")

files_limb = [f for f in os.listdir('../results/SVM_output/limb/') if f.find(".tsv") != -1 and f.find("Limb")!=-1 and f.find("nommus")==-1 and f.find("without") == -1 and f.find("highweightkmers") == -1]
files_nogc_limb = ["../results/SVM_output/limb/"+f for f in files_limb if f.find("gc") == -1 ]
files_gc_limb = ["../results/SVM_output/limb/"+f for f in files_limb if f.find("gc") != -1 ]

cmd_nogc = " ".join(["./preds2roc_pr_curves_customed.py","t"] + files_nogc_limb)
cmd_gc = " ".join(["./preds2roc_pr_curves_customed.py","t"] + files_gc_limb)

getaucs(cmd_nogc,limb=True)
getaucs(cmd_gc,gc="_gc",limb=True)

