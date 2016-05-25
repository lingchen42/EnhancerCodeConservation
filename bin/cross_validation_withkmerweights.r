#!/usr/bin/env Rscript
#Function: Train a species model(spectrum Kernel, k=5), perform 10 cross validation and output the scores. Default class size ration is 10.
#Usage:./cross_validation.r <species_name> <species_seqfile> <pos_num>

args <- commandArgs(trailingOnly = TRUE)
s <- args[1]
s_fh <- args[2]
s_posN <- as.integer(args[3])

library(caret)
library(kebabs)

dna <- readDNAStringSet(s_fh) 
labels <- rep.int(1, times=s_posN) #number of human enhancers
labels_neg <- rep.int(-1, times=length(dna)-s_posN) # number of negative regions
labels_all <- c(labels, labels_neg)
numSamples <- sample(1:length(dna))
flds <- createFolds(numSamples, k=10)

wts <- 10/c(1, 10)
names(wts) <- c(1,-1)
specK5 <- spectrumKernel(k=5, revComplement=FALSE)

for ( n in 1:10) {

test <- flds[[n]]
train <- c(1:length(dna))[-test]

model <- kbsvm(x=dna[train], y=labels_all[train], kernel=specK5, pkg="LibLineaR", svm="C-svc", cost=15, perfParameters="ALL", classWeights=wts) 
#extract feature weights of the model
fw <- t(sort(as.data.frame(featureWeights(model))))
out_name <- sprintf("%s_model_highweightkmers_run%s_%s.tsv",s,n,Sys.Date())
cat("Writing the feature weights to ", out_name)
write.table(fw, file=out_name,sep="\t")


  if ( n == 1 ) {
  names <- names(dna[test])
  true_labels <- labels_all[test]
  response <- predict(model, dna[test], raw=True, predictionType="response")
  decision <- predict(model, dna[test], raw=True, predictionType="decision") 
  fold <- rep.int(n, times=length(test))
  } else {
  names <- c(names, names(dna[test]))
  true_labels <- c(true_labels, labels_all[test])
  response <- c(response, predict(model, dna[test], raw=True, predictionType="response"))
  temp <- predict(model, dna[test], raw=True, predictionType="decision")
  decision <- c(decision, temp)
  rocdata <- computeROCandAUC(temp, labels_all[test], unique(labels_all[test]))
  fold <- c(fold, rep.int(n, times=length(test)))
  }

}

out_df <- data.frame(names, true_labels,response, decision, fold)
out_name <- sprintf("%s_cv_scores_%s.tsv", s, Sys.Date())
cat("Writing cross validation scores to ", out_name)
write.table(out_df, file=out_name, sep="\t", row.names=FALSE)




