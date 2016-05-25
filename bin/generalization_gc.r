#!/usr/bin/env Rscript
# Function: Use the model (spectrum kernel, k=5) trained on the first species seqs, apply to the second species. 
#           Output the model feature weights, the performance of first species model on second species. Default class size ratio is 10.
# Usage: ./generalization.r <species1_name> <species1_seqfile>  <pos_num> <species2_name> <species2_seqfile> <pos_num>

#parse arguments
args <- commandArgs(trailingOnly = TRUE)
s1 <- args[1]
s1_fh <- args[2]
s1_posN <- as.integer(args[3])
s2 <-args[4]
s2_fh <- args[5]
s2_posN <- as.integer(args[6])

library(kebabs)

species1 <- readDNAStringSet(s1_fh) #include all species1 enhancers and 10 random negative regions
labels_all <- c(rep.int(1, times=s1_posN), rep.int(-1, times=length(species1)-s1_posN) )
numSamples <- length(species1)

trainingFraction <- 1
train <- sample(1:numSamples, trainingFraction*numSamples)
test <- c(1:numSamples)[-train]
specK5 <- spectrumKernel(k=5, revComplement=FALSE) #spectrum kernel k=5

#class weights, 1:10
wts <- 10/c(1, 10) 
names(wts) <- c(1,-1)

#train a model based on species1 sequences
model <- kbsvm(x=species1[train], y=labels_all[train], kernel=specK5, pkg="LibLineaR", svm="C-svc", cost=15, perfParameters="ALL", classWeights=wts) 
#extract feature weights of the model
fw <- t(sort(as.data.frame(featureWeights(model))))
out_name <- sprintf("../results/SVM_output/%s_model_highweightkmers_%s_gc.tsv",s1,Sys.Date())
cat("Writing the feature weights to ", out_name)
write.table(fw, file=out_name,sep="\t")

#apply the model to species2, test on half of the data due to the size of data
species2 <- readDNAStringSet(s2_fh )
species2_labels <- c(rep.int(1, times=s2_posN ), rep.int(-1, times=length(species2)-s2_posN))
species2_test_fraction <- 0.5
species2_test1 <- sample(1:length(species2), species2_test_fraction*length(species2))

pred1 <- predict(model, species2[species2_test1])
preddec1 <- predict(model, species2[species2_test1], predictionType="decision")
evaluatePrediction(pred1,species2_labels[species2_test1], allLabels=unique(species2_labels[species2_test1]),decValues=preddec1)
names <- names(species2[species2_test1])
true_labels <- species2_labels[species2_test1]
response <- predict(model, species2[species2_test1], raw=True, predictionType="response")
decision <- predict(model, species2[species2_test1], raw=True, predictionType="decision")
out_df <- data.frame("names"= names, "labels"=true_labels, "response"=response, "decision"=decision)
out_name <- sprintf("../results/SVM_output/%s_applyto_%s_%s_gc.tsv", s1,s2,Sys.Date())
cat("Writing the prediction file of species1 trained model on species2 to ", out_name)
write.table(out_df, file=out_name, sep="\t", row.names=FALSE)
