#!usr/bin/env bash

#roadmap chromHMM annotated enhancers
mkdir ../results/overlapwithEn/
./getstates.py ../data/roadmap/roadmap_marks_inalltissues.bed ../results/overlapwithEn/roadmapEn.bed 6 7 12

sort -k1,1 -k2,2n ../results/overlapwithEn/roadmapEn.bed > ../results/overlapwithEn/sorted.bed
bedtools merge -i  ../results/overlapwithEn/sorted.bed > ../results/overlapwithEn/roadmapEn_merged.bed

rm  ../results/overlapwithEn/sorted.bed

#get FP TN FN TPs
# NOT GC CONTROLLED
./extractFPTN.py ../results/SVM_output/Hsap_cv_scores_2015-08-07.tsv 
#intersect with roadmap enhancers
wc -l ../results/overlapwithEn/Hsap_cv_scores_2015-08-07_FPs.bed
bedtools intersect -a ../results/overlapwithEn/Hsap_cv_scores_2015-08-07_FPs.bed -b ../results/overlapwithEn/roadmapEn_merged.bed -wa | sort -u | wc -l 
#68090 out of 73737 (92.34%)
wc -l ../results/overlapwithEn/Hsap_cv_scores_2015-08-07_TNs.bed
bedtools intersect -a ../results/overlapwithEn/Hsap_cv_scores_2015-08-07_TNs.bed -b ../results/overlapwithEn/roadmapEn_merged.bed -wa | sort -u | wc -l
#149982 out of 213941 (70.10%)
wc -l ../results/overlapwithEn/Hsap_cv_scores_2015-08-07_TPs.bed
bedtools intersect -a ../results/overlapwithEn/Hsap_cv_scores_2015-08-07_TPs.bed -b ../results/overlapwithEn/roadmapEn_merged.bed -wa | sort -u | wc -l 
#21683 out of 21703(99.91%)
wc -l ../results/overlapwithEn/Hsap_cv_scores_2015-08-07_FNs.bed
bedtools intersect -a ../results/overlapwithEn/Hsap_cv_scores_2015-08-07_FNs.bed -b ../results/overlapwithEn/roadmapEn_merged.bed -wa | sort -u | wc -l 
#7416/7449(99.56%)


#GC CONTROLLED
./extractFPTN.py ../results/SVM_output//Hsap_cv_scores_2015-09-22_gc.tsv
#intersect with roadmap enhancers
wc -l ../results/overlapwithEn/Hsap_cv_scores_2015-09-22_gc_FPs.bed
bedtools intersect -a ../results/overlapwithEn/Hsap_cv_scores_2015-09-22_gc_FPs.bed -b ../results/overlapwithEn/roadmapEn_merged.bed -wa | sort -u | wc -l 
#86575 out of 95083 (91.05%)
wc -l ../results/overlapwithEn/Hsap_cv_scores_2015-09-22_gc_TNs.bed
bedtools intersect -a ../results/overlapwithEn/Hsap_cv_scores_2015-09-22_gc_TNs.bed -b ../results/overlapwithEn/roadmapEn_merged.bed -wa | sort -u | wc -l 
#153898 out of 196437 (78.34%)
wc -l ../results/overlapwithEn/Hsap_cv_scores_2015-09-22_gc_TPs.bed
bedtools intersect -a ../results/overlapwithEn/Hsap_cv_scores_2015-09-22_gc_TPs.bed -b ../results/overlapwithEn/roadmapEn_merged.bed -wa | sort -u | wc -l 
#20526/20540.0 (99.93%)
wc -l ../results/overlapwithEn/Hsap_cv_scores_2015-09-22_gc_FNs.bed
bedtools intersect -a ../results/overlapwithEn/Hsap_cv_scores_2015-09-22_gc_FNs.bed -b ../results/overlapwithEn/roadmapEn_merged.bed -wa | sort -u | wc -l 
#8573/8612.0 (99.55%)
