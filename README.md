#Enhancer Code Conservation Documentation

We developed a machine-learning based cross-species enhancer prediction framwork to explore the conservation of enhancer sequence code in mammalian species.

Here we go through the analyses we performed in our manuscript. 

### Work Pipeline
<img>
### Data
We used two dataset, liver enhancers of 6 mammals (Villar et al. 2015) and limb enhancers of 3 mammals (Cotney et al. 2013). We then generated non-GC-controlled negatives and GC-controlled negatives for each set of enhancers. 
These bed regions can be found in the `data/seq/`

### Figure 2: 5-mer specturm SVM models can accurately identified enhancers in diverse mammals.(Within species enhancer prediction)

First, we demonstrated that the 5-mer specutrm SVM can accurately identified enhancers in diverse mammals.
We use `twoBitToFa` to convert the bed files to fasta sequences and use  `cross_validation.r` to perform within species enhancer predictions. The usage is
```
./cross_validation.r <Species> <Enhancers and Negatives sequence file> <Number of enhancers>
```

The complete commands we used to perform within species enhancer prediction are listed in `withinspecies.sh`.

We evaluated the within species enhancer prediction performace with ROC and PR curves.
```
#Figure 2a liver 
./preds2roc_pr_curves_customed.py  ../results/figures/liver_withinspecies_2a ../results/SVM_output/Hsap_cv_scores_2015-08-07.tsv ../results/SVM_output/Mmul_cv_scores_2015-08-07.tsv ../results/SVM_output/Mmus_cv_scores_2015-08-10.tsv ../results/SVM_output/Btau_cv_scores_2015-08-07.tsv ../results/SVM_output/Cfam_cv_scores_2015-08-07.tsv ../results/SVM_output/Mdom_cv_scores_2015-08-14.tsv
#Figure 2b limnb 
./preds2roc_pr_curves_customed.py  ../results/figures/limb_withinspecies_2b ../results/SVM_output/limb/Limb_Hsap_cv_scores_2016-01-07.tsv ../results/SVM_output/limb/Limb_Mmus_cv_scores_2016-01-07.tsv ../results/SVM_output/limb/Limb_Mmul_cv_scores_2016-01-07.tsv
```
### Figure 3: Enhancer classifiers perform well when applied across species (Cross species enhancer prediction)

Next, we applied SVM models trained in one species across species with `generalization.r`. The usage is
./generalization.r <species_1, used to trained the model> <species_1 Enhancers and Negatives sequence file> <Number of enhancers>  <species_2, used to trained the model> <species_2 Enhancers and Negatives sequence file> <Number of enhancers> 
The complete commands we used to perform cross species enhancer prediction are listed in `crossspecies.sh`.

The SVM models trained with human liver enhancers can predict enhancers in other species as well.
```
#Figure 3a
./preds2roc_pr_curves_customed.py ../results/figures/liver_human_applytoothers_3a ../results/SVM_output/Hsap_cv_scores_2015-08-07.tsv ../results/SVM_output/Hsap_applyto_Mmul_2015-11-19.tsv ../results/SVM_output/Hsap_applyto_Mmus_2015-11-19.tsv ../results/SVM_output/Hsap_applyto_Btau_2015-11-18.tsv ../results/SVM_output/Hsap_applyto_Cfam_2015-11-18.tsv ../results/SVM_output/Hsap_applyto_Mdom_2015-11-18.tsv
```
Same for the limb.
```
#Figure 3c
./preds2roc_pr_curves_customed.py ../results/figures/limb_human_applyttoothers_3c  ../results/SVM_output/limb/Limb_Hsap_cv_scores_2016-01-07.tsv ../results/SVM_output/limb/Limb_Hsap_applyto_Mmul_2016-01-06.tsv ../results/SVM_output/limb/Limb_Hsap_applyto_Mmus_2016-01-06.tsv
```
And, the good generalization of liver enhancer SVM models is true for every species pair, evaluated by computing auROC, auPR and relative auROC, relative auPR.
```
#Figure 3b and 3d
./get_auc.py
for FILE in ../results/figures/heatmaps/*tab
do
  ./heat_map_auc.py $FILE
done
```

### Figure 4: The predictions of enhancer classifiers trained in different species were strongly correlated.
We found that the predictions of enhancer classifiers trained in different species were strongly correlated. The commands for this analysis is in `./colorcode.sh`

In addition, we found that the mouse and opossum enhancers had a different GC content distribution, which contributed to the failure of mouse and opossum SVM model to recognize high GC content human enhancers.
```
./dist_gc.py  ../data/seqs/Hsap_enhancer_and_negative_regions/Hsap_Enhancers_gapexcluded_gcandlength.tsv ../data/seqs/Mmul_enhancer_and_negative_regions/Mmul_Enhancers_gapexcluded_gcandlength.tsv  ../data/seqs/Mmus_enhancer_and_negative_regions/Mmus_Enhancers_gapexcluded_gcandlength.tsv  ../data/seqs/Btau_enhancer_and_negative_regions/Btau_Enhancers_gapexcluded_gcandlength.tsv ../data/seqs/Cfam_enhancer_and_negative_regions/Cfam_Enhancers_gapexcluded_gcandlength.tsv ../data/seqs/Mdom_enhancer_and_negative_regions/Mdom_Enhancers_gapexcluded_gcandlength.tsv
```

### Figure 5: The ability to predict enhancers across species was maintained when the contribution from GC content is removed.

We then performed the GC-contolled analysis. The complete commands are in `crossspecies.sh`. And the performance was evaluated by computing auROC, auPR and relative auROC, relative auPR.
```
#Figure 5a and 5b
./get_auc.py
for FILE in ../results/figures/heatmaps/*gc*tab
do
  ./heat_map_auc.py $FILE
done
```

### Figure 6: The DNA sequence patterns most predictive of liver activity across species matched a common set of transcription factors.

We mapped the top 5-mers from liver enhancer species models of all species and found them match a common ser of transcription factors. The complete commands are in `Tfcompare.sh`

The Venn diagram is generated by [jVenn](http://bioinfo.genotoul.fr/jvenn/)
