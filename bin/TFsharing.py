#!/usr/bin/env python

import sys
import pandas as pd
import numpy as np

fns = sys.argv[1:]
species = []
tf_dict = dict()

for fn in fns:
  #make dictionary dt_dict: species:"top_pos":[], "neg_pos":[] ...
  if fn.find("Limb") == -1:
    s = fn.split("/")[-1].split("_")[0]
  else:
    s = fn.split("/")[-1].split("_")[1]
  species.append(s)

  df = pd.read_table(fn)
  cols = df.columns

  t_dict = dict()
  for c in cols:
    t_dict[c] = list(df[c].dropna()) 

  tf_dict[s] = t_dict

###############################################################################
#Compute shared in all
ALL = dict()

for i in species:
  for c in cols:
    try:
      ALL[c] = ALL[c].intersection(set(tf_dict[i][c]))
    except KeyError:
      ALL[c] = set(tf_dict[i][c])

all_ran = []
all_ranrc = []

for key in ALL.keys():
  if key.startswith("random"):
    if key.find("rc") ==-1:
      all_ran.append(len(ALL[key]))
    else:
      all_ranrc.append(len(ALL[key]))          
  elif key.find("rc") == -1:        
    print "ALL_"+key.upper()+"\t"+str(len(ALL[key]))
    print "ALL_"+key.upper()+"\t"+str(ALL[key])

print "ALL_RAN" +"\t"+str(np.mean(all_ran))+"+/-"+str(np.std(all_ran))
#print "ALL_RANRC" +"\t"+str(np.mean(all_ranrc))+"+/-"+str(np.std(all_ranrc))

###############################################################################
#PAIR WISE

# s1 s2
# names   s1  s2  shared
# top_pos
# top_neg 
## top_posrc
## top_negrc 
# rand  
## rand_rc

for index, i in enumerate(species):

  for j in species[index+1:]:

    print "\t".join(["","%s"%(i),"%s"%(j),"shared"])

    #compute random sharing of TFs from random set
    ran_i = []
    ran_j = []
    ran_count = []
    ran_i_rc = []
    ran_j_rc = []
    ran_count_rc = []

    for c in cols:

      if c.startswith("random"):
        if c.find("rc") ==-1:
          ran_i.append(len(set(tf_dict[i][c])))
          ran_j.append(len(set(tf_dict[j][c])))
          ran_count.append(len(set(tf_dict[i][c])&set(tf_dict[j][c])))
        else:
          ran_i_rc.append(len(set(tf_dict[i][c])))
          ran_j_rc.append(len(set(tf_dict[j][c])))
          ran_count_rc.append(len(set(tf_dict[i][c])&set(tf_dict[j][c])))

      elif c.find("rc") == -1:
        s1 =  set(tf_dict[i][c])
        s2 =  set(tf_dict[j][c])
        share = len(s1&s2)
        #print c.upper(),":",s1&s2
        print "\t".join([c.upper(),str(len(s1)),str(len(s2)),str(share)])

    print "\t".join(["RAN",str(np.mean(ran_i))+"+/-"+str(np.std(ran_i)),str(np.mean(ran_j))+"+/-"+str(np.std(ran_j)),str(np.mean(ran_count))+"+/-"+str(np.std(ran_count))])
    #print "\t".join(["RAN_RC",str(np.mean(ran_i_rc))"+/-"+str(np.std(ran_i_rc)),str(np.mean(ran_j_rc))+"+/-"+str(np.std(ran_j_rc)),str(np.mean(ran_count_rc))+"+/-"+str(np.std(ran_count_rc))])


