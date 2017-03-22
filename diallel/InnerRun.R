#!/usr/local/bin/Rscript

##---------------------------------------------------------------------------------------------------------------------
## Title: InnerRunner Script
## Author: Paul L. Maurizio
## Institution: UNC-Chapel Hill
## Program: Bioinformatics and Computational Biology Ph.D. Curriculum
## Advisors: Will Valdar, Mark Heise
## Date Created: 2016-04-06
## Date Updated: 2017-01-01
##---------------------------------------------------------------------------------------------------------------------

library(cmdline)
library(configfile)
library(tools)
library(methods)
library(treatmentResponseDiallel)

args <- c(cmdline.string("configfile"), cmdline.string("f"))
confile   	<- args[[1]]
config <- read.configfile(confile)
source("loadConfig.R")
datafile	<- args[[2]]
savedir			<- file_path_sans_ext(datafile)
dir.create(savedir, recursive=TRUE, showWarnings=FALSE)

sink(file.path(savedir, "output_innerrunner.out"))

##---------------------------------
## Run BayesDiallel
##---------------------------------

ZeroOutBeforeSample <- 0

if(batch=="NULL"){batch <- NULL}
if(fixed=="NULL"){fixed <- NULL}
if(random=="NULL"){random <- NULL}

random <- c(random,batch)
if(length(random)>1){random <- as.list(random)}
if(length(fixed)>1){fixed <- as.list(fixed)}

object	<- run.tr.diallel(filename=datafile, savedir=savedir,
			treatment=subset,
			trt.colname=trt_colname, 
			phenotype=phenotype, random=random, fixed=fixed, type=type, 
			BS=BS, burnin=burnin, thin=thin, lengthChains=lengthChains, 
			numChains=numChains)

sink()

cat("\n================================= \n")
cat("BayesDiallel complete for file=", basename(datafile), " \n", sep="")
cat(date(), "\n")
cat("================================= \n \n")

#sink(file.path(savedir, "output_innerrunner.out"), append=TRUE)
#inner.plotter(TreatDi=object[[2]], plotdir=savedir)
#sink()

if(!file.exists("InnerRun_Completed.txt")){
	try(system("touch InnerRun_Completed.txt", intern=TRUE))
}else{
	try(system("rm InnerRun_Completed.txt", intern=TRUE))
}
cat(datafile, "\n", file="InnerRun_Completed.txt", append=TRUE, sep="")

cat("\n================================= \n")
cat("InnerRun complete for file=", basename(datafile), " \n", sep="")
cat(date(), "\n")
cat("================================= \n \n")

system(paste("find ",  savedir, " -type f -name 'TreatDi*.RData' -delete &", sep=""))
system(paste("find ",  savedir, " -type f -name 'AutoDiallelChain*.bin' -delete &", sep=""))
system(paste("find ",  savedir, " -type f -name 'SaveTBSR5.RData' -delete &", sep=""))
system(paste("find ",  savedir, " -type f -name 'PartSaveListTBSR5.RData' -delete &", sep=""))
