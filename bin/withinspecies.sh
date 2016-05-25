#liver not gc controlled
#bed2fa
twoBitToFa -bed=../data/seqs/Hsap_enhancer_and_negative_regions/Hsap_enhancer_and_negative_region.bed -bedPos /dors/capra_lab/data/dna/hg19.2bit ../data/seqs/Hsap_enhancer_and_negative_regions/Hsap_enhancer_and_negative_region.fa&
twoBitToFa -bed=../data/seqs/Btau_enhancer_and_negative_regions/Btau_enhancer_and_negative_region.bed -bedPos /dors/capra_lab/data/dna/bosTau6.2bit ../data/seqs/Btau_enhancer_and_negative_regions/Btau_enhancer_and_negative_region.fa&
twoBitToFa -bed=../data/seqs/Mmus_enhancer_and_negative_regions/Mmus_enhancer_and_negative_region.bed -bedPos /dors/capra_lab/data/dna/mm10.2bit ../data/seqs/Mmus_enhancer_and_negative_regions/Mmus_enhancer_and_negative_region.fa&
twoBitToFa -bed=../data/seqs/Cfam_enhancer_and_negative_regions/Cfam_enhancer_and_negative_region.bed -bedPos /dors/capra_lab/data/dna/canFam3.2bit ../data/seqs/Cfam_enhancer_and_negative_regions/Cfam_enhancer_and_negative_region.fa&
twoBitToFa -bed=../data/seqs/Mmul_enhancer_and_negative_regions/Mmul_enhancer_and_negative_region.bed -bedPos /dors/capra_lab/data/dna/rheMac2/rheMac2.2bit ../data/seqs/Mmul_enhancer_and_negative_regions/Mmul_enhancer_and_negative_region.fa&
twoBitToFa -bed=../data/seqs/Mdom_enhancer_and_negative_regions/Mdom_enhancer_and_negative_region.bed -bedPos /dors/capra_lab/data/dna/monDom5/monDom5.2bit ../data/seqs/Mdom_enhancer_and_negative_regions/Mdom_enhancer_and_negative_region.fa&

#training
./cross_validation.r Btau ../data/seqs/Btau_enhancer_and_negative_regions/Btau_enhancer_and_negative_region.fa 30892
./cross_validation.r Mmus ../data/seqs/Mmus_enhancer_and_negative_regions/Mmus_enhancer_and_negative_region.fa 18517
./cross_validation.r Cfam ../data/seqs/Cfam_enhancer_and_negative_regions/Cfam_enhancer_and_negative_region.fa 18966
./cross_validation.r Mmul ../data/seqs/Mmul_enhancer_and_negative_regions/Mmul_enhancer_and_negative_region.fa 22911
./cross_validation.r Hsap ../data/seqs/Hsap_enhancer_and_negative_regions/Hsap_enhancer_and_negative_region.fa 29152
./cross_validation.r Mdom ../data/seqs/Mdom_enhancer_and_negative_regions/Mdom_enhancer_and_negative_region.fa 23160

#liver gc controlled
#bed2fa
twoBitToFa -bed=../data/seqs/Hsap_enhancer_and_negative_regions/Hsap_enhancer_and_negative_region_gc-aware.bed4 -bedPos /dors/capra_lab/data/dna/hg19.2bit ../data/seqs/Hsap_enhancer_and_negative_regions/Hsap_enhancer_and_negative_region_gc-aware.fa&
twoBitToFa -bed=../data/seqs/Btau_enhancer_and_negative_regions/Btau_enhancer_and_negative_region_gc-aware.bed4 -bedPos /dors/capra_lab/data/dna/bosTau6.2bit ../data/seqs/Btau_enhancer_and_negative_regions/Btau_enhancer_and_negative_region_gc-aware.fa&
twoBitToFa -bed=../data/seqs/Mmus_enhancer_and_negative_regions/Mmus_enhancer_and_negative_region_gc-aware.bed4 -bedPos /dors/capra_lab/data/dna/mm10.2bit ../data/seqs/Mmus_enhancer_and_negative_regions/Mmus_enhancer_and_negative_region_gc-aware.fa&
twoBitToFa -bed=../data/seqs/Cfam_enhancer_and_negative_regions/Cfam_enhancer_and_negative_region_gc-aware.bed4 -bedPos /dors/capra_lab/data/dna/canFam3.2bit ../data/seqs/Cfam_enhancer_and_negative_regions/Cfam_enhancer_and_negative_region_gc-aware.fa&
twoBitToFa -bed=../data/seqs/Mmul_enhancer_and_negative_regions/Mmul_enhancer_and_negative_region_gc-aware.bed4 -bedPos /dors/capra_lab/data/dna/rheMac2/rheMac2.2bit ../data/seqs/Mmul_enhancer_and_negative_regions/Mmul_enhancer_and_negative_region_gc-aware.fa&
twoBitToFa -bed=../data/seqs/Mdom_enhancer_and_negative_regions/Mdom_enhancer_and_negative_region_gc-aware.bed4 -bedPos /dors/capra_lab/data/dna/monDom5/monDom5.2bit ../data/seqs/Mdom_enhancer_and_negative_regions/Mdom_enhancer_and_negative_region_gc-aware.fa&
#training
./cross_validation.r Btau ../data/seqs/Btau_enhancer_and_negative_regions/Btau_enhancer_and_negative_region_gc-aware.fa 30892
./cross_validation.r Mmus ../data/seqs/Mmus_enhancer_and_negative_regions/Mmus_enhancer_and_negative_region_gc-aware.fa 18517
./cross_validation.r Cfam ../data/seqs/Cfam_enhancer_and_negative_regions/Cfam_enhancer_and_negative_region_gc-aware.fa 18966
./cross_validation.r Mmul ../data/seqs/Mmul_enhancer_and_negative_regions/Mmul_enhancer_and_negative_region_gc-aware.fa 22911
./cross_validation.r Hsap ../data/seqs/Hsap_enhancer_and_negative_regions/Hsap_enhancer_and_negative_region_gc-aware.fa 29152
./cross_validation.r Mdom ../data/seqs/Mdom_enhancer_and_negative_regions/Mdom_enhancer_and_negative_region_gc-aware.fa 23160

#limb not gc controlled
#bed2fa
twoBitToFa -bed=../data/seqs/Hsap_enhancer_and_negative_regions/Limb_Hsap_enhancer_and_negative_region_NEW.bed -bedPos /dors/capra_lab/data/dna/hg19.2bit ../data/seqs/Hsap_enhancer_and_negative_regions/Limb_Hsap_enhancer_and_negative_region.fa&
twoBitToFa -bed=../data/seqs/Mmus_enhancer_and_negative_regions/Limb_Mmus_enhancer_and_negative_region_NEW.bed -bedPos /dors/capra_lab/data/dna/mm9/mm9.2bit ../data/seqs/Mmus_enhancer_and_negative_regions/Limb_Mmus_enhancer_and_negative_region.fa&
twoBitToFa -bed=../data/seqs/Mmul_enhancer_and_negative_regions/Limb_Mmul_enhancer_and_negative_region_NEW.bed -bedPos /dors/capra_lab/data/dna/rheMac2/rheMac2.2bit ../data/seqs/Mmul_enhancer_and_negative_regions/Limb_Mmul_enhancer_and_negative_region.fa&

#training
./cross_validation.r Limb_Hsap ../data/seqs/Hsap_enhancer_and_negative_regions/Limb_Hsap_enhancer_and_negative_region.fa 25034
./cross_validation.r Limb_Mmul ../data/seqs/Mmul_enhancer_and_negative_regions/Limb_Mmul_enhancer_and_negative_region.fa 88560
./cross_validation.r Limb_Mmus ../data/seqs/Mmus_enhancer_and_negative_regions/Limb_Mmus_enhancer_and_negative_region.fa 87046

mv ../results/SVM_output/Limb* ../results/SVM_output/limb/


#limb gc controlled
#bed2fa
twoBitToFa -bed=../data/seqs/Hsap_enhancer_and_negative_regions/LimbHsap_enhancer_and_negative_region_gc-aware.bed -bedPos /dors/capra_lab/data/dna/hg19.2bit ../data/seqs/Hsap_enhancer_and_negative_regions/Limb_Hsap_enhancer_and_negative_region_gc-aware.fa&
twoBitToFa -bed=../data/seqs/Mmus_enhancer_and_negative_regions/LimbMmus_enhancer_and_negative_region_gc-aware.bed -bedPos /dors/capra_lab/data/dna/mm9/mm9.2bit ../data/seqs/Mmus_enhancer_and_negative_regions/Limb_Mmus_enhancer_and_negative_region_gc-aware.fa&
twoBitToFa -bed=../data/seqs/Mmul_enhancer_and_negative_regions/LimbMmul_enhancer_and_negative_region_gc-aware.bed -bedPos /dors/capra_lab/data/dna/rheMac2/rheMac2.2bit ../data/seqs/Mmul_enhancer_and_negative_regions/Limb_Mmul_enhancer_and_negative_region_gc-aware.fa&

#training
./cross_validation.r Limb_Hsap ../data/seqs/Hsap_enhancer_and_negative_regions/Limb_Hsap_enhancer_and_negative_region_gc-aware.fa 25034
./cross_validation.r Limb_Mmul ../data/seqs/Mmul_enhancer_and_negative_regions/Limb_Mmul_enhancer_and_negative_region_gc-aware.fa 88560
./cross_validation.r Limb_Mmus ../data/seqs/Mmus_enhancer_and_negative_regions/Limb_Mmus_enhancer_and_negative_region_gc-aware.fa 88560 
