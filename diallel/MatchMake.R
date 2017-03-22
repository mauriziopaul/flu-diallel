#!/usr/local/bin/Rscript

##---------------------------------------------------------------------------------------------------------------------
## Title: MatchMaker Script to Run BayesDiallel on Treatment Response Diallel Data
## Author: Paul L. Maurizio
## Institution: UNC-Chapel Hill
## Program: Bioinformatics and Computational Biology Ph.D. Curriculum
## Advisors: Will Valdar, Mark Heise
## Date Created: 2016-04-06
## Date Updated: 2017-01-01
##---------------------------------------------------------------------------------------------------------------------

library(cmdline)
library(configfile)
library(treatmentResponseDiallel)

cat("inside MatchMake\n")
args <- c(cmdline.string("configfile"), cmdline.string("datafile"))
cat("loaded args\n")
confile   	<- args[[1]]
datafile	<- args[[2]]
config <- read.configfile(confile)
cat("read configfile\n")
source("loadConfig.R")

##---------------------------------
## LOAD DATA, IMPUTE MATCHED PAIRS
##---------------------------------

data <- read.csv(datafile)

if(!(subset=="All")){
		data <- data[data[subset_colname]==subset,]
	}

reptotal <- paste0("NumMatches_", sprintf("%03d",reps))
fdir <- savedir

timethis(filenames	<- make.matches(data=data, reps=reps, trt.string=trt_string, 
	ctrl.string=ctrl_string, fdir=fdir, strategy=strategy, select.phenotype=phenotype,
	force=FALSE))

cat(filenames, sep="\n", file="MIMP_List.txt", append=FALSE)

if(!file.exists(file.path(fdir, "PostRunConfig.sh"))){
	try(system("touch PostRunConfig.sh", intern = TRUE))
	try(system(paste("touch ", fdir, "/PostRunConfig.sh", sep=""), intern = TRUE))
	cat("DIR=\'", fdir, "\'\n", file=file.path(fdir,"PostRunConfig.sh"), append=FALSE, sep="")
	cat("THIN=", thin, "\n", file=file.path(fdir,"PostRunConfig.sh"), append=TRUE, sep="")

	cat("DIR=\'", fdir, "\'\n", file="PostRunConfig.sh", append=FALSE, sep="")
	cat("THIN=", thin, "\n", file="PostRunConfig.sh", append=TRUE, sep="")
}

cat("\n================================= \n")
cat("Matches made. (n=", reps, ") \n", sep="")
cat("PostRunConfig.sh saved.")
cat(date(), "\n")
cat("================================= \n \n")

