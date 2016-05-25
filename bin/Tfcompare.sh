#!/usr/bin/env bash

#Liver GC controlled
#generate kmer list
./kmer_list.py  ../results/SVM_output/Hsap_model_highweightkmers_2015-11-19_gc.tsv
./kmer_list.py  ../results/SVM_output/Btau_model_highweightkmers_2015-11-18_gc.tsv
./kmer_list.py ../results/SVM_output/Cfam_model_highweightkmers_2015-11-19_gc.tsv
./kmer_list.py ../results/SVM_output/Mmul_model_highweightkmers_2015-11-19_gc.tsv
./kmer_list.py  ../results/SVM_output/Mmus_model_highweightkmers_2015-11-18_gc.tsv
./kmer_list.py ../results/SVM_output/Mdom_model_highweightkmers_2015-11-19_gc.tsv

#covert kmer list to meme format
for FILE in ../data/kmer/kebabs_kmer_list/*2015-11*gc*list
do
./kmer2meme_general.py $FILE
done

#JASPAR CORE 2014 VERTEBRATES
for FILE in ../data/kmer/meme_format_kmers/*model*2015-11*gc**.meme
do
FN=$(basename $FILE)
mkdir -p ../results/tomtom/JASPAR_CORE_2014_vertebrates/$FN/
tomtom -query-pseudo 0.001 -oc ../results/tomtom/JASPAR_CORE_2014_vertebrates/$FN/ $FILE ../data/motifdatabase/JASPAR_CORE_2014_vertebrates.meme
done

#Hsap
./kmermatchedTFtable.py ../results/tomtom/JASPAR_CORE_2014_vertebrates/Hsap_model_highweightkmers_2015-11-19_gc.meme/tomtom.txt ../results/SVM_output/Hsap_model_highweightkmers_2015-11-19_gc.tsv  ../data/TFexpression/JASPAR_CORE_2014_vertebrates_expression_profile.txt ../results/TFanalysis/matched_TFs_liver_GCcontrolled/
#Mmul
./kmermatchedTFtable.py ../results/tomtom/JASPAR_CORE_2014_vertebrates/Mmul_model_highweightkmers_2015-11-19_gc.meme/tomtom.txt ../results/SVM_output/Mmul_model_highweightkmers_2015-11-19_gc.tsv ../data/TFexpression/JASPAR_CORE_2014_vertebrates_expression_profile.txt ../results/TFanalysis/matched_TFs_liver_GCcontrolled/
#Mmus
./kmermatchedTFtable.py ../results/tomtom/JASPAR_CORE_2014_vertebrates/Mmus_model_highweightkmers_2015-11-18_gc.meme/tomtom.txt ../results/SVM_output/Mmus_model_highweightkmers_2015-11-18_gc.tsv ../data/TFexpression/JASPAR_CORE_2014_vertebrates_expression_profile.txt ../results/TFanalysis/matched_TFs_liver_GCcontrolled/
#Btau
./kmermatchedTFtable.py ../results/tomtom/JASPAR_CORE_2014_vertebrates/Btau_model_highweightkmers_2015-11-18_gc.meme/tomtom.txt ../results/SVM_output/Btau_model_highweightkmers_2015-11-18_gc.tsv  ../data/TFexpression/JASPAR_CORE_2014_vertebrates_expression_profile.txt ../results/TFanalysis/matched_TFs_liver_GCcontrolled/
#Cfam
./kmermatchedTFtable.py ../results/tomtom/JASPAR_CORE_2014_vertebrates/Cfam_model_highweightkmers_2015-11-19_gc.meme/tomtom.txt ../results/SVM_output/Cfam_model_highweightkmers_2015-11-19_gc.tsv ../data/TFexpression/JASPAR_CORE_2014_vertebrates_expression_profile.txt ../results/TFanalysis/matched_TFs_liver_GCcontrolled/
#Mdom
./kmermatchedTFtable.py ../results/tomtom/JASPAR_CORE_2014_vertebrates/Mdom_model_highweightkmers_2015-11-19_gc.meme/tomtom.txt ../results/SVM_output/Mdom_model_highweightkmers_2015-11-19_gc.tsv ../data/TFexpression/JASPAR_CORE_2014_vertebrates_expression_profile.txt ../results/TFanalysis/matched_TFs_liver_GCcontrolled/

./TFsharing.py ../results/TFanalysis/matched_TFs_liver_GCcontrolled/Hsap_model_highweightkmers_2015-11-19_gc_TFmatches.tab ../results/TFanalysis/matched_TFs_liver_GCcontrolled/Btau_model_highweightkmers_2015-11-18_gc_TFmatches.tab ../results/TFanalysis/matched_TFs_liver_GCcontrolled/Cfam_model_highweightkmers_2015-11-19_gc_TFmatches.tab ../results/TFanalysis/matched_TFs_liver_GCcontrolled/Mmul_model_highweightkmers_2015-11-19_gc_TFmatches.tab ../results/TFanalysis/matched_TFs_liver_GCcontrolled/Mmus_model_highweightkmers_2015-11-18_gc_TFmatches.tab ../results/TFanalysis/matched_TFs_liver_GCcontrolled/Mdom_model_highweightkmers_2015-11-19_gc_TFmatches.tab> ../results/TFanalysis/matched_TFs_liver_GCcontrolled/sharing.txt

###############################################################################################################################################################################################################
#Liver not GC controlled

#generate kmer list
./kmer_list.py  ../results/SVM_output/Hsap_model_highweightkmers_2015-11-19.tsv
./kmer_list.py   ../results/SVM_output/Btau_model_highweightkmers_2015-11-19.tsv
./kmer_list.py ../results/SVM_output/Cfam_model_highweightkmers_2015-11-19.tsv
./kmer_list.py ../results/SVM_output/Mmul_model_highweightkmers_2015-11-19.tsv
./kmer_list.py  ../results/SVM_output/Mmus_model_highweightkmers_2015-11-19.tsv
./kmer_list.py ../results/SVM_output/Mdom_model_highweightkmers_2015-11-19.tsv

#covert kmer list to meme format
for FILE in ../data/kmer/kebabs_kmer_list/*2015-11-19.list
do
./kmer2meme_general.py $FILE
done

#JASPAR CORE 2014 VERTEBRATES
for FILE in ../data/kmer/meme_format_kmers/*2015-11-19.meme
do
FN=$(basename $FILE)
mkdir -p ../results/tomtom/JASPAR_CORE_2014_vertebrates/$FN/
tomtom -query-pseudo 0.001 -oc ../results/tomtom/JASPAR_CORE_2014_vertebrates/$FN/  $FILE ../data/motifdatabase/JASPAR_CORE_2014_vertebrates.meme
done


#Hsap
./kmermatchedTFtable.py ../results/tomtom/JASPAR_CORE_2014_vertebrates/Hsap_model_highweightkmers_2015-11-19.meme/tomtom.txt ../results/SVM_output/Hsap_model_highweightkmers_2015-11-19.tsv ../data/TFexpression/JASPAR_CORE_2014_vertebrates_expression_profile.txt ../results/TFanalysis/matched_TFs_liver_notGC/
#Mmul
./kmermatchedTFtable.py ../results/tomtom/JASPAR_CORE_2014_vertebrates/Mmul_model_highweightkmers_2015-11-19.meme/tomtom.txt ../results/SVM_output/Mmul_model_highweightkmers_2015-11-19.tsv ../data/TFexpression/JASPAR_CORE_2014_vertebrates_expression_profile.txt ../results/TFanalysis/matched_TFs_liver_notGC/
#Mmus
./kmermatchedTFtable.py ../results/tomtom/JASPAR_CORE_2014_vertebrates/Mmus_model_highweightkmers_2015-11-19.meme/tomtom.txt ../results/SVM_output/Mmus_model_highweightkmers_2015-11-19.tsv ../data/TFexpression/JASPAR_CORE_2014_vertebrates_expression_profile.txt ../results/TFanalysis/matched_TFs_liver_notGC/
#Btau
./kmermatchedTFtable.py ../results/tomtom/JASPAR_CORE_2014_vertebrates/Btau_model_highweightkmers_2015-11-19.meme/tomtom.txt ../results/SVM_output/Btau_model_highweightkmers_2015-11-19.tsv ../data/TFexpression/JASPAR_CORE_2014_vertebrates_expression_profile.txt ../results/TFanalysis/matched_TFs_liver_notGC/
#Cfam
./kmermatchedTFtable.py ../results/tomtom/JASPAR_CORE_2014_vertebrates/Cfam_model_highweightkmers_2015-11-19.meme/tomtom.txt ../results/SVM_output/Cfam_model_highweightkmers_2015-11-19.tsv ../data/TFexpression/JASPAR_CORE_2014_vertebrates_expression_profile.txt ../results/TFanalysis/matched_TFs_liver_notGC/
#Mdom
./kmermatchedTFtable.py ../results/tomtom/JASPAR_CORE_2014_vertebrates/Mdom_model_highweightkmers_2015-11-19.meme/tomtom.txt ../results/SVM_output/Mdom_model_highweightkmers_2015-11-19.tsv ../data/TFexpression/JASPAR_CORE_2014_vertebrates_expression_profile.txt ../results/TFanalysis/matched_TFs_liver_notGC/

./TFsharing.py ../results/TFanalysis/matched_TFs_liver_notGC/Hsap_model_highweightkmers_2015-11-19_TFmatches.tab ../results/TFanalysis/matched_TFs_liver_notGC/Btau_model_highweightkmers_2015-11-19_TFmatches.tab ../results/TFanalysis/matched_TFs_liver_notGC/Cfam_model_highweightkmers_2015-11-19_TFmatches.tab ../results/TFanalysis/matched_TFs_liver_notGC/Mmul_model_highweightkmers_2015-11-19_TFmatches.tab ../results/TFanalysis/matched_TFs_liver_notGC/Mmus_model_highweightkmers_2015-11-19_TFmatches.tab  ../results/TFanalysis/matched_TFs_liver_notGC/Mdom_model_highweightkmers_2015-11-19_TFmatches.tab> ../results/TFanalysis/matched_TFs_liver_notGC/sharing.txt

###############################################################################################################################################################################################################
#Limb not GC controlled

#generate kmer list
./kmer_list.py ../results/SVM_output/limb/Limb_Hsap_model_highweightkmers_2016-01-06.tsv
./kmer_list.py ../results/SVM_output/limb/Limb_Mmul_model_highweightkmers_2016-01-06.tsv
./kmer_list.py ../results/SVM_output/limb/Limb_Mmus_model_highweightkmers_2016-01-06.tsv

#covert kmer list to meme format
for FILE in ../data/kmer/kebabs_kmer_list/Limb*list
do
./kmer2meme_general.py $FILE
done

#JASPAR CORE 2014 VERTEBRATES
for FILE in ../data/kmer/meme_format_kmers/Limb*meme
do
FN=$(basename $FILE)
mkdir -p ../results/tomtom/JASPAR_CORE_2014_vertebrates/$FN/
tomtom -query-pseudo 0.001 -oc ../results/tomtom/JASPAR_CORE_2014_vertebrates/$FN/  $FILE ../data/motifdatabase/JASPAR_CORE_2014_vertebrates.meme
done

#Hsap
./kmermatchedTFtable.py ../results/tomtom/JASPAR_CORE_2014_vertebrates/Limb_Hsap_model_highweightkmers_2016-01-06.meme/tomtom.txt ../results/SVM_output/limb/Limb_Hsap_model_highweightkmers_2016-01-06.tsv ../data/TFexpression/JASPAR_CORE_2014_vertebrates_expression_profile.txt ../results/TFanalysis/matched_TFs_limb_notGC/
#Mmul
./kmermatchedTFtable.py ../results/tomtom/JASPAR_CORE_2014_vertebrates/Limb_Mmul_model_highweightkmers_2016-01-06.meme/tomtom.txt ../results/SVM_output/limb/Limb_Mmul_model_highweightkmers_2016-01-06.tsv ../data/TFexpression/JASPAR_CORE_2014_vertebrates_expression_profile.txt ../results/TFanalysis/matched_TFs_limb_notGC/
#Mmus
./kmermatchedTFtable.py ../results/tomtom/JASPAR_CORE_2014_vertebrates/Limb_Mmus_model_highweightkmers_2016-01-06.meme/tomtom.txt ../results/SVM_output/limb/Limb_Mmus_model_highweightkmers_2016-01-06.tsv ../data/TFexpression/JASPAR_CORE_2014_vertebrates_expression_profile.txt ../results/TFanalysis/matched_TFs_limb_notGC/

./TFsharing.py ../results/TFanalysis/matched_TFs_limb_notGC/Limb_Hsap_model_highweightkmers_2016-01-06_TFmatches.tab ../results/TFanalysis/matched_TFs_limb_notGC/Limb_Mmul_model_highweightkmers_2016-01-06_TFmatches.tab ../results/TFanalysis/matched_TFs_limb_notGC/Limb_Mmus_model_highweightkmers_2016-01-06_TFmatches.tab  > ../results/TFanalysis/matched_TFs_limb_notGC/sharing.txt

###############################################################################################################################################################################################################
#Limb GC controlled
#generate kmer list
./kmer_list.py  ../results/SVM_output/limb/Limb_Hsap_model_highweightkmers_2016-01-18_gc.tsv
./kmer_list.py ../results/SVM_output/limb/Limb_Mmul_model_highweightkmers_2016-01-18_gc.tsv
./kmer_list.py ../results/SVM_output/limb/Limb_Mmus_model_highweightkmers_2016-01-18_gc.tsv

#covert kmer list to meme format
for FILE in ../data/kmer/kebabs_kmer_list/Limb*gc*list
do
./kmer2meme_general.py $FILE
done

#JASPAR CORE 2014 VERTEBRATES
for FILE in ../data/kmer/meme_format_kmers/Limb*gc*meme
do
FN=$(basename $FILE)
mkdir -p ../results/tomtom/JASPAR_CORE_2014_vertebrates/$FN/
tomtom -query-pseudo 0.001 -oc ../results/tomtom/JASPAR_CORE_2014_vertebrates/$FN/  $FILE ../data/motifdatabase/JASPAR_CORE_2014_vertebrates.meme
done

#Hsap
./kmermatchedTFtable.py ../results/tomtom/JASPAR_CORE_2014_vertebrates/Limb_Hsap_model_highweightkmers_2016-01-18_gc.meme/tomtom.txt ../results/SVM_output/limb/Limb_Hsap_model_highweightkmers_2016-01-18_gc.tsv ../data/TFexpression/JASPAR_CORE_2014_vertebrates_expression_profile.txt ../results/TFanalysis/matched_TFs_limb_GCcontrolled/
#Mmul
./kmermatchedTFtable.py ../results/tomtom/JASPAR_CORE_2014_vertebrates/Limb_Mmul_model_highweightkmers_2016-01-18_gc.meme/tomtom.txt ../results/SVM_output/limb/Limb_Mmul_model_highweightkmers_2016-01-18_gc.tsv ../data/TFexpression/JASPAR_CORE_2014_vertebrates_expression_profile.txt ../results/TFanalysis/matched_TFs_limb_GCcontrolled/
#Mmus
./kmermatchedTFtable.py ../results/tomtom/JASPAR_CORE_2014_vertebrates/Limb_Mmus_model_highweightkmers_2016-01-18_gc.meme/tomtom.txt ../results/SVM_output/limb/Limb_Mmus_model_highweightkmers_2016-01-18_gc.tsv ../data/TFexpression/JASPAR_CORE_2014_vertebrates_expression_profile.txt ../results/TFanalysis/matched_TFs_limb_GCcontrolled/

./TFsharing.py ../results/TFanalysis/matched_TFs_limb_GCcontrolled/Limb_Hsap_model_highweightkmers_2016-01-18_gc_TFmatches.tab ../results/TFanalysis/matched_TFs_limb_GCcontrolled/Limb_Mmul_model_highweightkmers_2016-01-18_gc_TFmatches.tab ../results/TFanalysis/matched_TFs_limb_GCcontrolled/Limb_Mmus_model_highweightkmers_2016-01-18_gc_TFmatches.tab  > ../results/TFanalysis/matched_TFs_limb_GCcontrolled/sharing.txt

