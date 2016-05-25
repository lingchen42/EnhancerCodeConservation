#!/usr/bin/env python

import pandas as pd
from scipy import stats
import sys
import numpy as np

df1 = pd.read_table(sys.argv[1],header=None)
df2 = pd.read_table(sys.argv[2],header=None)

print "The lower, the broader of activity"
print "FILE1 mean: ",np.mean(df1[1])
print "FILE2 mean: ",np.mean(df2[1])
#s,p=stats.ttest_ind(df1[1],df2[1], equal_var = False)
s,p=stats.mannwhitneyu(df1[1],df2[1])

#print "From unequal var t-test: statistic ",s," p-value:", p 
print "From mannwhitneyu test: statistic ",s," p-value:", p
