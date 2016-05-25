#!/usr/bin/env python
#usage kmer2meme <kmer list/a file with a list of kmers>

import sys
import os
import string
import random

def kmer2matrix(kmer):
  matrix=[['0' for x in range(len(kmer))] for x in range(4)] 
  for index, i in enumerate(kmer):
    if i == "A":
      matrix[0][index]='1'
    if i == "C":
      matrix[1][index]='1'
    if i == "G":
      matrix[2][index]='1'
    if i == "T":
      matrix[3][index]='1'
  return matrix

def writekmer(kmer):  
  matrix=kmer2matrix(kmer)        
  with open(''.join([temp_dir, kmer,'.pfm']),'w+') as outfh:
    for x in matrix:
      outfh.write('\t'.join(x)+'\n')

OUTDIR = "../data/kmer/meme_format_kmers/"
if not os.path.exists(OUTDIR): os.makedirs(OUTDIR)

temp_dir="temp_%s/"%(''.join(random.choice(string.digits) for _ in range(8)))
if not os.path.exists(temp_dir): 
  os.makedirs(temp_dir)

try: 
  fh = open(sys.argv[1],"r")
  
  kmers = []
  for line in fh.readlines():
    kmers.append(line[:-1])

  out = OUTDIR+sys.argv[1].split("/")[-1].replace("list","meme")
  print "generate ", OUTDIR+sys.argv[1].split("/")[-1].replace("list","meme")
except IOError:
  kmers = sys.argv[1:]
  out = OUTDIR+"_".join(kmers)+".meme"
  print "generate ", OUTDIR+"_".join(kmers)+".meme"


for kmer in kmers:
  writekmer(kmer)

cmd1="jaspar2meme -pfm %s > %s" %(temp_dir, out)

os.system(cmd1)
cmd2="rm -r %s"%(temp_dir)
os.system(cmd2)


