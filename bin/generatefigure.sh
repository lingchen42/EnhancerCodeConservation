#!/usr/bin/env bash

## FIGURE 2
#within species training
#sbatch withinspecies.slurm

#data visualization
#2a liver withinspecies
./preds2roc_pr_curves_customed.py  ../results/figures/liver_withinspecies_1a ../results/SVM_output/Hsap_cv_scores_2015-08-07.tsv ../results/SVM_output/Mmul_cv_scores_2015-08-07.tsv ../results/SVM_output/Mmus_cv_scores_2015-08-10.tsv ../results/SVM_output/Btau_cv_scores_2015-08-07.tsv ../results/SVM_output/Cfam_cv_scores_2015-08-07.tsv ../results/SVM_output/Mdom_cv_scores_2015-08-14.tsv
#2b limnb withinsepcies
./preds2roc_pr_curves_customed.py  ../results/figures/limb_withinspecies_1b ../results/SVM_output/limb/Limb_Hsap_cv_scores_2016-01-07.tsv ../results/SVM_output/limb/Limb_Mmus_cv_scores_2016-01-07.tsv ../results/SVM_output/limb/Limb_Mmul_cv_scores_2016-01-07.tsv

####################################################################################################################################################################################
## FIGURE 2
#cross species training
#sbatch crossspecies.slurm

#data visualization
# #2a Human liver classifier can do well across species 
./preds2roc_pr_curves_customed.py ../results/figures/liver_human_applytoothers_2a ../results/SVM_output/Hsap_cv_scores_2015-08-07.tsv ../results/SVM_output/Hsap_applyto_Mmul_2015-11-19.tsv ../results/SVM_output/Hsap_applyto_Mmus_2015-11-19.tsv ../results/SVM_output/Hsap_applyto_Btau_2015-11-18.tsv ../results/SVM_output/Hsap_applyto_Cfam_2015-11-18.tsv ../results/SVM_output/Hsap_applyto_Mdom_2015-11-18.tsv

# #2b,4,s3,s4 liver,limb pair-wise analysis 
./get_auc.py
for FILE in ../results/figures/heatmaps/*tab
do
  ./heat_map_auc.py $FILE
done

## s2 sup removal of conserved enhancers does not affect the performance
./preds2roc_pr_curves_customed.py ../results/figures/noconserveden_s1 ../results/SVM_output/Hsap_nommus_applyto_Mmus_nohsap_2015-11-19.tsv ../results/SVM_output/nommus_Hsap_cv_scores_2015-11-20.tsv

# #2c Human limb classifer  can do well across species
./preds2roc_pr_curves_customed.py ../results/figures/limb_human_applyttoothers_2c  ../results/SVM_output/limb/Limb_Hsap_cv_scores_2016-01-07.tsv ../results/SVM_output/limb/Limb_Hsap_applyto_Mmul_2016-01-06.tsv ../results/SVM_output/limb/Limb_Hsap_applyto_Mmus_2016-01-06.tsv

####################################################################################################################################################################################
#Figure 3 & s
./colorcode.sh&
#Figure 3f GC  content distribution
./dist_gc.py  ../data/seqs/Hsap_enhancer_and_negative_regions/Hsap_Enhancers_gapexcluded_gcandlength.tsv ../data/seqs/Mmul_enhancer_and_negative_regions/Mmul_Enhancers_gapexcluded_gcandlength.tsv  ../data/seqs/Mmus_enhancer_and_negative_regions/Mmus_Enhancers_gapexcluded_gcandlength.tsv  ../data/seqs/Btau_enhancer_and_negative_regions/Btau_Enhancers_gapexcluded_gcandlength.tsv ../data/seqs/Cfam_enhancer_and_negative_regions/Cfam_Enhancers_gapexcluded_gcandlength.tsv ../data/seqs/Mdom_enhancer_and_negative_regions/Mdom_Enhancers_gapexcluded_gcandlength.tsv

####################################################################################################################################################################################
#Figure 4 generates along with figure 2b

####################################################################################################################################################################################
#Figure 5

#kmer sharing
#liver not gc controlled 5%, no rc
mkdir ../results/figures/venn/topkmerlist/
./compare_kmers_norc.py 0.05  ../results/SVM_output/Hsap_model_highweightkmers_2015-11-19.tsv  ../results/SVM_output/Mmul_model_highweightkmers_2015-11-19.tsv  ../results/SVM_output/Mmus_model_highweightkmers_2015-11-19.tsv  ../results/SVM_output/Btau_model_highweightkmers_2015-11-19.tsv  ../results/SVM_output/Cfam_model_highweightkmers_2015-11-19.tsv  ../results/SVM_output/Mdom_model_highweightkmers_2015-11-19.tsv
#liver gc controlled 5%,no rc
./compare_kmers_norc.py 0.05  ../results/SVM_output/Hsap_model_highweightkmers_2015-11-19_gc.tsv  ../results/SVM_output/Mmul_model_highweightkmers_2015-11-19_gc.tsv  ../results/SVM_output/Mmus_model_highweightkmers_2015-11-18_gc.tsv  ../results/SVM_output/Btau_model_highweightkmers_2015-11-18_gc.tsv  ../results/SVM_output/Cfam_model_highweightkmers_2015-11-19_gc.tsv  ../results/SVM_output/Mdom_model_highweightkmers_2015-11-19_gc.tsv
head -n 51  ../results/figures/venn/topkmerlist/*gc*pos*  >  ../results/figures/venn/topposkmers_liver_gc.txt
#limb not gc contolled
mkdir ../results/figures/venn/topkmerlist_limb/
./compare_kmers_norc.py 0.05 ../results/SVM_output/limb/Limb_Hsap_model_highweightkmers_2016-01-06.tsv ../results/SVM_output/limb/Limb_Mmus_model_highweightkmers_2016-01-06.tsv ../results/SVM_output/limb/Limb_Mmul_model_highweightkmers_2016-01-06.tsv
head -n 51  ../results/figures/venn/topkmerlist_limb/*pos*  >  ../results/figures/venn/topposkmers_limb.txt

#limb gc contolled
./compare_kmers_norc.py 0.05 ../results/SVM_output/limb/Limb_Hsap_model_highweightkmers_2016-01-18_gc.tsv ../results/SVM_output/limb/Limb_Mmus_model_highweightkmers_2016-01-18_gc.tsv ../results/SVM_output/limb/Limb_Mmul_model_highweightkmers_2016-01-18_gc.tsv
head -n 51  ../results/figures/venn/topkmerlist_limb/*gc*pos*  >  ../results/figures/venn/topposkmers_limb_gc.txt

#TF sharing
./Tfcompare.sh

#Go to http://bioinfo.genotoul.fr/jvenn/example.html, generate the graph and download the data in ../results/figures/venn/
