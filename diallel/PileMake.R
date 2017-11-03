#!/usr/local/bin/Rscript

##---------------------------------------------------------------------------------------------------------------------
## Title: PileMaker Script
## Author: Paul L. Maurizio
## Institution: UNC-Chapel Hill
## Program: Bioinformatics and Computational Biology Ph.D. Curriculum
## Advisors: Will Valdar, Mark Heise
## Date Created: 2016-04-06
## Date Updated: 2017-01-01
##---------------------------------------------------------------------------------------------------------------------

library(cmdline)
library(configfile)

args <- c(cmdline.string("configfile"), cmdline.string("MIMPlist"))
confile   	<- args[[1]]
MIMPlist	<- args[[2]]
config <- read.configfile(confile)
source("loadConfig.R")

MIMPfiles <- scan(file.path(savedir,MIMPlist), what="", sep="\n")

cat("#!/bin/bash \n", file="BatchSend.sh", append=FALSE)

pilename_list <- NULL

for(p in 1:num_piles){
	this_pilename <- paste("Pile", sprintf("%02d",p),".sh", sep="")
	cat("#!/bin/bash \n", file=this_pilename, append=FALSE)
	cat("bash ", this_pilename, " & \n", sep="", file="BatchSend.sh", append=TRUE)
	pilename_list <- c(pilename_list, this_pilename)
}

for(i in c(1:length(MIMPfiles))){
	j <- i%%num_piles + 1
	this_pilename <- as.character(pilename_list[j])
	this_mimpfilename <- as.character(MIMPfiles[i])
	cat("Rscript InnerRun.R --args --configfile=", confile, " --f=\'", this_mimpfilename, 
		"\' && \n", sep="", file=this_pilename, append=TRUE)
}

for(p in 1:num_piles){
	cat("\n================================= \n")
	cat("echo \"Finished pile ", pilename_list[p], "\" \n", sep="", file=pilename_list[p], append=TRUE)
	cat(date(), "\n")
	cat("================================= \n \n")
}

cat("\n================================= \n")
cat("Piles made. (n=", num_piles, ") \n", sep="")
cat(date(), "\n")
cat("================================= \n \n")

