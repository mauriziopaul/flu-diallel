#!/usr/local/bin/Rscript

##---------------------------------------------------------------------------------------------------------------------
## Title: Match Maker Script for Matching Treated and Control Animals in the Diallel
## Author: Paul L. Maurizio
## Institution: UNC-Chapel Hill
## Program: Bioinformatics and Computational Biology Ph.D. Curriculum
## Advisors: Will Valdar, Mark Heise
## Date Created: 2016-05-29
## Date Updated: 2017-01-01
##---------------------------------------------------------------------------------------------------------------------

library(cmdline)
library(configfile)
library(treatmentResponseDiallel)
data(PiximusData)

cat("Inside ImputedMatch\n")
confile <- cmdline.string("configfile")
cat("Loaded args\n")
config <- read.configfile(confile)
cat("Configfile read\n")
source("loadConfig.R")
covariate <- covImpute

# cat(filenames, sep="\n", file=file.path(savedir, "Imputed_List.txt"), append=FALSE)
filenames <- as.vector(read.delim(file=file.path(savedir, "Imputed_List.txt"), sep="\n", header=FALSE)[[1]])

matchnames <- NULL
for(i in 1:reps){
  matchname <- file.path(savedir, paste("TreatDi_MP", sprintf("%04d",i), ".csv", sep=""))
  this.matchname <- make.matches(data=read.csv(filenames[i]), reps=1, trt.string=trt_string, 
    ctrl.string=ctrl_string, fdir=savedir, strategy=strategy, select.phenotype=phenotype,
    force=FALSE, matchname=matchname)
  matchnames <- c(as.character(matchnames), as.character(this.matchname))
  cat("Imputed Data Set #:" , i, "\t")
}

cat(as.character(matchnames), sep="\n", file=file.path(savedir, "MIMQ_List.txt"), append=FALSE)

if(!file.exists(file.path(savedir, "PostRunConfig.sh"))){
  try(system("touch PostRunConfig.sh", intern = TRUE))
  try(system(paste("touch ", savedir, "/PostRunConfig.sh", sep=""), intern = TRUE))
  cat("DIR=\'", savedir, "\'\n", file=file.path(savedir,"PostRunConfig.sh"), append=FALSE, sep="")
  cat("THIN=", thin, "\n", file=file.path(savedir,"PostRunConfig.sh"), append=TRUE, sep="")

  cat("DIR=\'", savedir, "\'\n", file="PostRunConfig.sh", append=FALSE, sep="")
  cat("THIN=", thin, "\n", file="PostRunConfig.sh", append=TRUE, sep="")
}


cat("\n================================= \n")
cat("Matches made. (n=", reps, ") \n", sep="")
cat("PostRunConfig.sh saved.")
cat(date(), "\n")
cat("================================= \n \n")


## Plot matched data to see how estimates change across imputations

pdf(file.path(savedir,"matched_data.pdf"), width=6, height=6)

for(i in 1:reps){
  dat.mat <- read.csv(matchnames[i])
  plot(dat.mat[,phenotype]~dat.mat[,covariate], data=dat.mat, col="black", pch=16, 
	ylim=range(dat.mat[,phenotype]), xlim=range(dat.mat[,covariate]))
  points(dat.mat[,phenotype]~dat.mat[,covariate], data=dat.mat, col="white", pch=16, cex=0.7)
  cat("Done plotting rep:", i, "\n")
}

dev.off()

