#!/usr/local/bin/Rscript

##---------------------------------------------------------------------------------------------------------------------
## Title: Posterior Prediction Script to Impute Missing Values for Treatment Response Diallel Data
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

cat("Inside PostPredImpute\n")
cat("Loaded args\n")
confile     <- cmdline.string("configfile")
config <- read.configfile(confile)
cat("Configfile read\n")
source("loadConfig.R")
covariate <- covImpute

#dir.create(savedir, recursive=TRUE)
system(paste("mkdir ", savedir, sep=""))
system(paste("cp ", confile, " ", savedir, "/", sep=""))

dat <- read.csv(datafile)

dat.sub <- dat[dat[trt_colname]==trt_string,]

if(phenbatch=="NULL"){phenbatch <- NULL}
if(phenfixed=="NULL"){phenfixed <- NULL}
if(phenrandom=="NULL"){phenrandom <- NULL}
if(covbatch=="NULL"){covbatch <- NULL}
if(covfixed=="NULL"){covfixed <- NULL}
if(covrandom=="NULL"){covrandom <- NULL}


##---------------------------
## Impute Missing Covariates
##---------------------------

## We will impute covariate weights for all "missing" animals
## using all the mice in the study.
## We assume the mice were selected for treatment at random, and therefore the covariates should
## be distributed the same between treatment groups.
## We will still model batch when imputing the covariates.

AFD.cov <- DiallelAnalyzer(data = dat, father.strain="sire",
                           mother.strain="dam", phenotype=covImpute, is.female="is_female",
                           FixedEffects=covfixed,
                           RandomEffects=c(covbatch, covrandom),
                           Models=Example.Piximus.Models[1], sigmasq.start=1,  
                           numChains=numChains, lengthChains=lengthChains,
                           burnin=1,
                           DIC.Only=FALSE, tauPriorFile=Example.Piximus.tau.Prior.Info,
                           SaveAFDFile=file.path(savedir, "CovPostPredImpute.RData"),
                           LogTransform=FALSE, DoBayesSpike=FALSE, DoFirstCenter=TRUE);

# We will need Sigma and Batch/Block effects, at each timestep
var.labels.all <- varnames(AFD.cov$AllDiallelObs[[1]]$cent.chains)
new.labels <- sapply(X=var.labels.all, FUN=var.translate, USE.NAMES=FALSE)
effectEst.cov <- mcmc.stack.and.burn(AFD.cov$AllDiallelObs[[1]]$cent.chains, burnin=500)
varnames(effectEst.cov) <- new.labels
batchnames.cov <- grep("^Batch:", colnames(effectEst.cov), value = TRUE)
batchnames.cov <- c(grep("^RandomEffect:", colnames(effectEst.cov), value = TRUE), batchnames.cov)
BlockEst.cov <- as.matrix(effectEst.cov[,batchnames.cov])
SigmaEst.cov <- as.vector(effectEst.cov[,"Sigma"])
#FixedEst.cov <- as.vector(effectEst.cov[,"FixedEffect:1"])

# We will also need posterior predictive means at each time step
ADO.cov <- AFD.cov$AllDiallelObs[[1]]
MyPostSum.cov <- PosteriorPredSummary(ADO.cov, AFD.cov, burnin = 1, AFD=AFD.cov, keep = TRUE);
postPred.cov <- mcmc.stack.and.burn(ADO.cov$PostKeeper$FakeCoda, burnin=500)

# Read in imputation table (dam, sire, sex, bat(ch))
toImpute <- read.csv(imputefile, stringsAsFactors=FALSE)
toImpute$batname <- paste("Batch:1:", toImpute$bat, sep="")
toImpute$Block <- toImpute$bat
toImpute$Strain <- paste(toImpute$dam, toImpute$sire, sep="")
toImpute$Day <- assignment
toImpute$Trt <- trt_string
toImpute[,"is_female"] <- ifelse(toImpute$Sex=="F", 1, ifelse(toImpute$Sex=="M", 0, NA))
toImpute$Mx1Dam   <- Mx1HaploClass(toImpute$dam)
toImpute$Mx1Sire  <- Mx1HaploClass(toImpute$sire)
toImpute$Mx1Diplo <- paste(toImpute$Mx1Dam, toImpute$Mx1Sire, sep="_")
toImpute$Mx1Diplo6 <- paste(apply(toImpute, MARGIN=1,
                        FUN=function(x){paste(sort(c(x["Mx1Dam"], x["Mx1Sire"])), collapse="_")}))
toImpute$catname <- paste(  "S:", toImpute[,"is_female"],
                           ";i:", lettersToNumbers(toImpute$dam),
                           ";k:", lettersToNumbers(toImpute$sire), sep="")

# Reduce toImpute to batches that are in data set, when dataset is subsetted
# toImpute <- droplevels(subset(toImpute), Block %in% dat$Block)

toCensor <- droplevels(subset(toImpute, action=="-"))
toImpute <- droplevels(subset(toImpute, action=="+"))

##-------------------------------
## Impute Missing Phenotype Data
##-------------------------------

AFD <- DiallelAnalyzer(data = dat.sub, father.strain="sire",
        mother.strain="dam", phenotype=phenImpute, is.female="is_female",
        FixedEffects=phenfixed,
        RandomEffects=c(phenbatch,phenrandom),
        Models=Example.Piximus.Models[1], sigmasq.start=1,  
        numChains=numChains, lengthChains=lengthChains,
        burnin=1,
        DIC.Only=FALSE, tauPriorFile=Example.Piximus.tau.Prior.Info,
        SaveAFDFile=file.path(savedir, "PhenPostPredImpute.RData"),
        LogTransform=FALSE, DoBayesSpike=FALSE, DoFirstCenter=TRUE);

# We will need Sigma and Batch/Block effects, as well as fixed effects, at each timestep
var.labels.all <- varnames(AFD$AllDiallelObs[[1]]$cent.chains)
new.labels <- sapply(X=var.labels.all, FUN=var.translate, USE.NAMES=FALSE)
effectEst <- mcmc.stack.and.burn(AFD$AllDiallelObs[[1]]$cent.chains, burnin=500)
varnames(effectEst) <- new.labels
batchnames <- grep("^Batch:", colnames(effectEst), value = TRUE)
batchnames <- c(grep("^RandomEffect:", colnames(effectEst), value = TRUE), batchnames)
BlockEst <- as.matrix(effectEst[,batchnames])
SigmaEst <- as.vector(effectEst[,"Sigma"])
FixedEst <- as.vector(effectEst[,"FixedEffect:1"])

# We will also need posterior predictive means at each time step
ADO <- AFD$AllDiallelObs[[1]]
MyPostSum <- PosteriorPredSummary(ADO, AFD, burnin = 1, AFD=AFD, keep = TRUE);
postPred <- mcmc.stack.and.burn(ADO$PostKeeper$FakeCoda, burnin=burnin)

##---------------
## ImputeMissing
##---------------

colsToCopy <- c("Block", "dam", "sire", "Strain", "Sex", "Day", "Trt", "is_female",
                "Mx1Dam", "Mx1Sire", "Mx1Diplo", "Mx1Diplo6") #, "D0"

timethis(filenames <- imputeMissing(dat=dat, phenotype=phenImpute, covariate=covImpute, 
              postPred=postPred, postPred.cov=postPred.cov,
              effectEst=effectEst, SigmaEst=SigmaEst, BlockEst=BlockEst, FixedEst=FixedEst,
              effectEst.cov=effectEst.cov, SigmaEst.cov=SigmaEst.cov, BlockEst.cov=BlockEst.cov,
              toImpute=toImpute, toCensor=toCensor, colsToCopy=colsToCopy,
              numImps=reps, savedir=savedir))

cat(as.character(filenames), sep="\n", file=file.path(savedir, "Imputed_List.txt"), append=FALSE)

### Stack imputed rows to see their distribution.
ncolumns <- length(names(dat))
imputed.check <- as.data.frame(matrix(data=rep(NA, times=ncolumns), nrow=1, ncol=ncolumns))
names(imputed.check) <- names(dat)
for(i in 1:reps){
  this.filename <- file.path(savedir, sprintf("%04d",i), "ImputedData.csv")
  thisdat <- read.csv(this.filename)
  imputed.check <- rbind(imputed.check, subset(thisdat, is.na(X)))
  cat("Finished rep: ", i,"\n")
}

dat.all <- read.csv(datafile)

pdf(file.path(savedir,"imputed_data.pdf"), width=6, height=6)

for(i in 1:reps){
  dat.imp <- read.csv(file.path(savedir, sprintf("%04d",i), "ImputedData.csv"))
  dat.sub <- subset(dat.imp, is.na(X))
  dat.subc <- subset(dat.imp, !is.na(X))
  
  plot(dat.all[,phenotype]~dat.all[,covariate], data=dat.all, col="black", pch=16,
    ylim=range(dat.all[,phenotype], na.rm=TRUE), xlim=range(dat.all[,covariate], na.rm=TRUE))
  points(dat.subc[,phenotype]~dat.subc[,covariate], data=dat.subc, col="white", pch=16, cex=0.7)
  points(dat.sub[,phenotype]~dat.sub[,covariate], data=dat.sub, col="black", pch=16)
  points(dat.sub[,phenotype]~dat.sub[,covariate], data=dat.sub, col="red", pch=16, cex=0.7)
  cat("Done plotting rep:", i, "\n")
}

dev.off()

write.csv(imputed.check, file.path(savedir, "ImputeCensorStack.csv"))

cat("\n================================= \n")
cat("Imputations made. (n=", reps, ") \n", sep="")
cat(date(), "\n")
cat("================================= \n \n")
