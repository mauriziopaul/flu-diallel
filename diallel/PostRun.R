#!/usr/local/bin/Rscript

##---------------------------------------------------------------------------------------------------------------------
## Title: PostRunner Script
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
library(tools)

confile <- cmdline.string("configfile")
config <- read.configfile(confile)
source("loadConfig.R")

MIMPfiles <- scan(file.path(savedir,MIMPlist), what="", sep="\n")
mimpdir <- file_path_sans_ext(MIMPfiles)

if(batch=="NULL"){batch <- NULL}
if(fixed=="NULL"){fixed <- NULL}
if(random=="NULL"){random <- NULL}

random <- c(random,batch)
if(length(random)>1){random <- as.list(random)}
if(length(fixed)>1){fixed <- as.list(fixed)}

sink(file.path(savedir, "output_postrunner.out"))

#stacked.matches <- as.mcmc(read.csv(file.path(savedir, "output/stackedTreatDiMIMPs.csv.gz"), 
stacked.matches <- as.mcmc(read.csv(file.path(savedir, "output/stackedTreatDiMIMPs.csv"), 
	check.names=FALSE))
#stacked.posterior <- read.delim2(file.path(savedir, "output/stackedPostTreatDiMIMPs.csv.gz"), 
stacked.posterior <- read.delim2(file.path(savedir, "output/stackedPostTreatDiMIMPs.csv"), 
	check.names=FALSE, sep=";", dec=",")
#stacked.postpred <- as.mcmc(read.csv(file.path(savedir, "output/stackedPostPredTreatDiMIMPs.csv.gz"), 
stacked.postpred <- as.mcmc(read.csv(file.path(savedir, "output/stackedPostPredTreatDiMIMPs.csv"), 
	check.names=FALSE))

##-----------------------------------------------
## Analyze Stacked
##-----------------------------------------------

## HPD PLOTS

timethis(
	hpd.plotter(plotdir=savedir, chain.object=stacked.matches, batched="TRUE", fixed="TRUE")
	)

## TRACE PLOTS

tracevars <- length(colnames(stacked.matches))
tracedir <- file.path(savedir, "output", "traceplots")
dir.create(tracedir, recursive=TRUE)

for(i in 1:tracevars){
tracename <- paste0("trace-", sprintf("%04d", i), 
	"-", colnames(stacked.matches)[i], ".pdf")
pdf(file.path(tracedir, tracename), width=10, height=5)
#png(file.path(tracedir, tracename), width = 1000, height = 500)
par(mfrow=c(1,2))
plot(as.vector(stacked.matches[,colnames(stacked.matches)[i]]), type="l", 
	main=colnames(stacked.matches)[i], ylab="")
plot(density(as.vector(stacked.matches[,colnames(stacked.matches)[i]])), type="l", 
	main=colnames(stacked.matches)[i], ylab="")
dev.off()
}

##------------
## PLOT VARPs
##------------
stacked.posterior <- as.mcmc(stacked.posterior)
VarPs <- summary(stacked.posterior)
VarPs1 <- VarPs[[1]]
VarPs2 <- VarPs[[2]]
colnames(VarPs2)[c(3,1,5)] <- c("mean", "lower.bound", "upper.bound") ## mean should be median
write.csv(VarPs2, file=file.path(savedir, "output/PSqTableStacked.csv"), row.names=TRUE)
try(plotVarps(fdir=savedir, fname="output/PSqTableStacked.csv")) ##

##------------
## PLOT MIPs
##------------
mip.filenames <- lapply(X=mimpdir, FUN=function(x){return(file.path(x, "MIP.csv"))})
try(mip.summary <- averageMips(filenames=mip.filenames, fdir=file.path(savedir,"output")))
try(plotMips(fdir=file.path(savedir, "output"), fname=mip.summary))

cat("NumChains: ", numChains, ", lengthChains: ", lengthChains,
	", Burnin: ", burnin, ", Thinning: ", thin, "\n", sep="")
cat(date(), "\n")

sink()

##----------
## PLOT POE
##----------

timethis( 
	## This step can be slow
	poe.analyzer(stacked.matches, plotdir=file.path(savedir,"output"))
)

cat("Finished plotting.\n")
