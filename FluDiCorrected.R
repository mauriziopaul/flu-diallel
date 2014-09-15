#!/usr/bin/Rscript

##---------------------------------------
## Title: Script to correct flu-infected estimates by mock estimates
## Author: Paul L. Maurizio
## Institution: UNC-Chapel Hill
## Program: Bioinformatics and Computational Biology Ph.D. Curriculum
## Advisors: Will Valdar, Mark Heise
## Date Created: 2014-08-21
## Date Updated: 2014-09-12
##
## To Do: (1) Plot according to PsyDiallel
##---------------------------------------

## rm(list=ls()); gc()
args <- commandArgs(TRUE)
good.arg <- NULL
good.arg[1] <- ".RData"
good.arg[2] <- ".Rdata"

try(require(stringr, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);

## Default setting when no arguments passed
if(!(length(args) > 6) | 
	(!((str_sub(args[1], start=-6) == good.arg[1]) | (str_sub(args[2], start=-6) == good.arg[2]))))
	{
		  stop("Please provide proper *.Rdata or *.RData file names")
}

load(file=args[1])
load(file=args[2])

print(commandArgs())
print(try(system("who", intern = TRUE)))
print(try(system("env", intern = TRUE)))

## library("lattice"); library("coda"); library("R.methodsS3"); 
## library("R.oo"); library("corpcor"); library("BayesDiallel")

try(require(lattice, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
try(require(coda, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
try(require(R.methodsS3, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
try(require(R.oo, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
try(require(corpcor, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE);
require(BayesDiallel, quietly = TRUE, warn.conflicts=FALSE)

data(PiximusData)

data <- read.csv("Flu_Diallel_v007.csv")
is_female <- sapply(data$Sex, function(x){ifelse(x=="F", 1, 0)})
##D4_diff <- data$D4-data$D0

## Correct data$D4 column for FLU animals

summary(AFD$AllDiallelObs[[1]]$centered.chains)[[1]]

pct_d4 <- 100*data$D4/data$D0
data <- cbind(data, is_female, pct_d4)
data.mock <- subset(data, Trt=="MOCK")
data.flu <- subset(data, Trt=="FLU")
data.mock.D4 <- subset(data, Trt=="FLU")
data.flu.D4 <- subset(data, Trt=="FLU")
print(head(data.flu))

## For Day 0, both mock and flu data can be analyzed together
## For percent starting weight, measure separate data sets
## For difference in starting weight, measure separate data sets

if(args[1]=="All"){
	data <- data
	}else{
		if(args[1]=="Flu"){
			data <- data.flu
			}else{
				data <- data.mock
			}
		}

head(data)
fname <- paste(args[1], args[2], args[3], args[4], sep="_")

FluDi = DiallelAnalyzer(data = data,
    father.strain="sire",
        mother.strain="dam",
        phenotype=args[2],
        is.female="is_female",
        FixedEffects=if(args[3]=="NULL"){NULL}else{"D0"},
        RandomEffects=if(args[4]=="NULL"){NULL}else{"Block"},
        Models=Example.Piximus.Models[1],
        sep="", na.strings="NA",
        sigmasq.start = 1, numChains =5, lengthChains = 3000,
        Verbose = TRUE,
        burnin = 500, thin = 1, DIC.Only = FALSE,
        tauPriorFile = Example.Piximus.tau.Prior.Info, SaveAFDFile=paste(fname, ".RData", sep=""), LogTransform=FALSE)

##save.image("FluDi.RData")
##save(AFD=FluDi, this=FluDi, file="FluDi.RData")
##load("FluDi.RData")

##names(FluDi)
##names(FluDi$AllDiallelObs[[1]])
##summary(FluDi$AllDiallelObs[[1]]$centered.chains)

##----------------------------------------
## ALAN TRANSLATOR, FROM WILL
##----------------------------------------
## source("AlanTranslator.R")

dir.create(fname)
pdf(file.path(fname,"HPD_plots.pdf"), width=11, height=8.5)
par(mfrow=c(1,4))
xlim=c(-10,10)
#plot(FluDi$AllDiallelObs[[length(FluDi$AllDiallelObs)]]$raw.chains)
#plot.hpd(FluDi$AllDiallelObs[[1]]$centered.chains)

## write(varnames(FluDi$AllDiallelObs[[1]]$raw.chains[,c(3,58:65, 74:81, 90:97)]),"")
## write()
var.labels <- c("BetaHybrid:Av", "aj:1 - mean(aj)", "aj:2 - mean(aj)", "aj:3 - mean(aj)", "aj:4 - mean(aj)", 
	      		"aj:5 - mean(aj)", "aj:6 - mean(aj)", "aj:7 - mean(aj)", "aj:8 - mean(aj)",
	      		"motherj:1 - mean(motherj)", "motherj:2 - mean(motherj)", "motherj:3 - mean(motherj)", "motherj:4 - mean(motherj)", 
				"motherj:5 - mean(motherj)", "motherj:6 - mean(motherj)", "motherj:7 - mean(motherj)", "motherj:8 - mean(motherj)",
				"dominancej:1 - mean(dominancej)", "dominancej:2 - mean(dominancej)", "dominancej:3 - mean(dominancej)", "dominancej:4 - mean(dominancej)",
				"dominancej:5 - mean(dominancej)", "dominancej:6 - mean(dominancej)", "dominancej:7 - mean(dominancej)", "dominancej:8 - mean(dominancej)")

new.labels <- c("inbreed.overall", "additive:AJ", "additive:B6", "additive:129", "additive:NOD", "additive:NZO", "additive:CAST",
       	  		"additive:PWK", "additive:WSB", 
			    "maternal:AJ", "maternal:B6", "maternal:129", "maternal:NOD", "maternal:NZO", "maternal:CAST",
			    "maternal:PWK", "maternal:WSB",
			    "inbreeding:AJ", "inbreeding:B6", "inbreeding:129", "inbreeding:NOD", "inbreeding:NZO", "inbreeding:CAST",
			    "inbreeding:PWK", "inbreeding:WSB")
plot.hpd(FluDi$AllDiallelObs[[1]]$centered.chains[,var.labels], names=new.labels, xlim=xlim, main="General effects") ## inbreeding, additive, maternal
abline(v=0, col="gray")

# write(varnames(FluDi$AllDiallelObs[[1]]$centered.chains[,c(106:133, 162:189)]), "")

var.labels <- c("SymCrossjk:j:2;k:1 - mean(SymCrossjk)", "SymCrossjk:j:3;k:1 - mean(SymCrossjk)", "SymCrossjk:j:4;k:1 - mean(SymCrossjk)", 
	      	"SymCrossjk:j:5;k:1 - mean(SymCrossjk)", "SymCrossjk:j:6;k:1 - mean(SymCrossjk)", "SymCrossjk:j:7;k:1 - mean(SymCrossjk)", 
		"SymCrossjk:j:8;k:1 - mean(SymCrossjk)", "SymCrossjk:j:3;k:2 - mean(SymCrossjk)", "SymCrossjk:j:4;k:2 - mean(SymCrossjk)", 
		"SymCrossjk:j:5;k:2 - mean(SymCrossjk)", "SymCrossjk:j:6;k:2 - mean(SymCrossjk)", "SymCrossjk:j:7;k:2 - mean(SymCrossjk)", 
		"SymCrossjk:j:8;k:2 - mean(SymCrossjk)", "SymCrossjk:j:4;k:3 - mean(SymCrossjk)", "SymCrossjk:j:5;k:3 - mean(SymCrossjk)", 
		"SymCrossjk:j:6;k:3 - mean(SymCrossjk)", "SymCrossjk:j:7;k:3 - mean(SymCrossjk)", "SymCrossjk:j:8;k:3 - mean(SymCrossjk)", 
		"SymCrossjk:j:5;k:4 - mean(SymCrossjk)", "SymCrossjk:j:6;k:4 - mean(SymCrossjk)", "SymCrossjk:j:7;k:4 - mean(SymCrossjk)", 
		"SymCrossjk:j:8;k:4 - mean(SymCrossjk)", "SymCrossjk:j:6;k:5 - mean(SymCrossjk)", "SymCrossjk:j:7;k:5 - mean(SymCrossjk)", 
		"SymCrossjk:j:8;k:5 - mean(SymCrossjk)", "SymCrossjk:j:7;k:6 - mean(SymCrossjk)", "SymCrossjk:j:8;k:6 - mean(SymCrossjk)", 
		"SymCrossjk:j:8;k:7 - mean(SymCrossjk)", 
		"ASymCrossjkDkj:j:2;k:1 - mean(ASymCrossjkDkj)", "ASymCrossjkDkj:j:3;k:1 - mean(ASymCrossjkDkj)", "ASymCrossjkDkj:j:4;k:1 - mean(ASymCrossjkDkj)", 
		"ASymCrossjkDkj:j:5;k:1 - mean(ASymCrossjkDkj)", "ASymCrossjkDkj:j:6;k:1 - mean(ASymCrossjkDkj)", "ASymCrossjkDkj:j:7;k:1 - mean(ASymCrossjkDkj)", 
		"ASymCrossjkDkj:j:8;k:1 - mean(ASymCrossjkDkj)", "ASymCrossjkDkj:j:3;k:2 - mean(ASymCrossjkDkj)", "ASymCrossjkDkj:j:4;k:2 - mean(ASymCrossjkDkj)", 
		"ASymCrossjkDkj:j:5;k:2 - mean(ASymCrossjkDkj)", "ASymCrossjkDkj:j:6;k:2 - mean(ASymCrossjkDkj)", "ASymCrossjkDkj:j:7;k:2 - mean(ASymCrossjkDkj)", 
		"ASymCrossjkDkj:j:8;k:2 - mean(ASymCrossjkDkj)", "ASymCrossjkDkj:j:4;k:3 - mean(ASymCrossjkDkj)", "ASymCrossjkDkj:j:5;k:3 - mean(ASymCrossjkDkj)", 
		"ASymCrossjkDkj:j:6;k:3 - mean(ASymCrossjkDkj)", "ASymCrossjkDkj:j:7;k:3 - mean(ASymCrossjkDkj)", "ASymCrossjkDkj:j:8;k:3 - mean(ASymCrossjkDkj)", 
		"ASymCrossjkDkj:j:5;k:4 - mean(ASymCrossjkDkj)", "ASymCrossjkDkj:j:6;k:4 - mean(ASymCrossjkDkj)", "ASymCrossjkDkj:j:7;k:4 - mean(ASymCrossjkDkj)", 
		"ASymCrossjkDkj:j:8;k:4 - mean(ASymCrossjkDkj)", "ASymCrossjkDkj:j:6;k:5 - mean(ASymCrossjkDkj)", "ASymCrossjkDkj:j:7;k:5 - mean(ASymCrossjkDkj)", 
		"ASymCrossjkDkj:j:8;k:5 - mean(ASymCrossjkDkj)", "ASymCrossjkDkj:j:7;k:6 - mean(ASymCrossjkDkj)", "ASymCrossjkDkj:j:8;k:6 - mean(ASymCrossjkDkj)", 
		"ASymCrossjkDkj:j:8;k:7 - mean(ASymCrossjkDkj)")

new.labels <- c("v:B6;AJ", "v:129;AJ", "v:NOD;AJ", "v:NZO;AJ", "v:CAST;AJ", "v:PWK;AJ", "v:WSB;AJ", 
	      		"v:129;B6", "v:NOD;B6", "v:NZO;B6", "v:CAST;B6", "v:PWK;B6", "v:WSB;B6", "v:NOD;129", 
				"v:NZO;129", "v:CAST;129", "v:PWK;129", "v:WSB;129", "v:NZO;NOD", "v:CAST;NOD", "v:PWK;NOD", 
				"v:WSB;NOD", "v:CAST;NZO", "v:PWK;NZO", "v:WSB;NZO", "v:PWK;CAST", "v:WSB;CAST", "v:WSB;PWK", 
				"w:B6;AJ", "w:129;AJ", "w:NOD;AJ", "w:NZO;AJ", "w:CAST;AJ", "w:PWK;AJ", "w:WSB;AJ", "w:129;B6", 
				"w:NOD;B6", "w:NZO;B6", "w:CAST;B6", "w:PWK;B6", "w:WSB;B6", "w:NOD;129", "w:NZO;129", "w:CAST;129", 
				"w:PWK;129", "w:WSB;129", "w:NZO;NOD", "w:CAST;NOD", "w:PWK;NOD", "w:WSB;NOD", "w:CAST;NZO", 
				"w:PWK;NZO", "w:WSB;NZO", "w:PWK;CAST", "w:WSB;CAST", "w:WSB;PWK")

plot.hpd(FluDi$AllDiallelObs[[1]]$centered.chains[,var.labels], names=new.labels, xlim=xlim, main="Strainpair-specific") ## symmetric, asymmetric
## plot.hpd(FluDi$AllDiallelObs[[1]]$centered.chains[,c(106:133, 162:189)]) ## symmetric, asymmetric
abline(v=0, col="gray")

## write(varnames(FluDi$AllDiallelObs[[1]]$centered.chains[,c(2,4,66:73, 82:89, 98:105)], "")
var.labels <- c( "Gender:Av", "BetaHybrid:Gender:Av", "Gender:aj:1 - mean(Gender:aj)", "Gender:aj:2 - mean(Gender:aj)", 
				"Gender:aj:3 - mean(Gender:aj)", "Gender:aj:4 - mean(Gender:aj)", "Gender:aj:5 - mean(Gender:aj)", 
				"Gender:aj:6 - mean(Gender:aj)", "Gender:aj:7 - mean(Gender:aj)", "Gender:aj:8 - mean(Gender:aj)", 
				"Gender:motherj:1 - mean(Gender:motherj)", "Gender:motherj:2 - mean(Gender:motherj)", "Gender:motherj:3 - mean(Gender:motherj)", 
				"Gender:motherj:4 - mean(Gender:motherj)", "Gender:motherj:5 - mean(Gender:motherj)", "Gender:motherj:6 - mean(Gender:motherj)", 
				"Gender:motherj:7 - mean(Gender:motherj)", "Gender:motherj:8 - mean(Gender:motherj)", "Gender:dominancej:1 - mean(Gender:dominancej)", 
				"Gender:dominancej:2 - mean(Gender:dominancej)", "Gender:dominancej:3 - mean(Gender:dominancej)", 
				"Gender:dominancej:4 - mean(Gender:dominancej)", "Gender:dominancej:5 - mean(Gender:dominancej)", 
				"Gender:dominancej:6 - mean(Gender:dominancej)", "Gender:dominancej:7 - mean(Gender:dominancej)", 
				"Gender:dominancej:8 - mean(Gender:dominancej)")
new.labels <- c("female.overall", "female.inbred", "additive:AJ", "additive:B6", "additive:129", "additive:NOD", "additive:NZO", "additive:CAST",
       	  		"additive:PWK", "additive:WSB",
       	  		"maternal:AJ", "maternal:B6", "maternal:129", "maternal:NOD", "maternal:NZO", "maternal:CAST",
			    "maternal:PWK", "maternal:WSB",
			    "inbreeding:AJ", "inbreeding:B6", "inbreeding:129", "inbreeding:NOD", "inbreeding:NZO", "inbreeding:CAST",
			    "inbreeding:PWK", "inbreeding:WSB")

plot.hpd(FluDi$AllDiallelObs[[1]]$centered.chains[,c(var.labels)], names=new.labels, xlim=xlim, main="Sex-specific") ## female overall, female inbred, female-specific add, mat, inbred
#plot.hpd(FluDi$AllDiallelObs[[1]]$centered.chains[,c(2, 4, 66:73, 82:89, 98:105)], names=new.labels) ## female overall, female inbred, female-specific add, mat, inbred
abline(v=0, col="gray")

## write(varnames(FluDi$AllDiallelObs[[1]]$centered.chains[,c(134:161, 190:217)]), "")
var.labels <- c("Gender:SymCrossjk:j:2;k:1 - mean(Gender:SymCrossjk)", "Gender:SymCrossjk:j:3;k:1 - mean(Gender:SymCrossjk)", 
				"Gender:SymCrossjk:j:4;k:1 - mean(Gender:SymCrossjk)", "Gender:SymCrossjk:j:5;k:1 - mean(Gender:SymCrossjk)", 
				"Gender:SymCrossjk:j:6;k:1 - mean(Gender:SymCrossjk)", "Gender:SymCrossjk:j:7;k:1 - mean(Gender:SymCrossjk)", 
				"Gender:SymCrossjk:j:8;k:1 - mean(Gender:SymCrossjk)", "Gender:SymCrossjk:j:3;k:2 - mean(Gender:SymCrossjk)", 
				"Gender:SymCrossjk:j:4;k:2 - mean(Gender:SymCrossjk)", "Gender:SymCrossjk:j:5;k:2 - mean(Gender:SymCrossjk)", 
				"Gender:SymCrossjk:j:6;k:2 - mean(Gender:SymCrossjk)", "Gender:SymCrossjk:j:7;k:2 - mean(Gender:SymCrossjk)", 
				"Gender:SymCrossjk:j:8;k:2 - mean(Gender:SymCrossjk)", "Gender:SymCrossjk:j:4;k:3 - mean(Gender:SymCrossjk)", 
				"Gender:SymCrossjk:j:5;k:3 - mean(Gender:SymCrossjk)", "Gender:SymCrossjk:j:6;k:3 - mean(Gender:SymCrossjk)", 
				"Gender:SymCrossjk:j:7;k:3 - mean(Gender:SymCrossjk)", "Gender:SymCrossjk:j:8;k:3 - mean(Gender:SymCrossjk)", 
				"Gender:SymCrossjk:j:5;k:4 - mean(Gender:SymCrossjk)", "Gender:SymCrossjk:j:6;k:4 - mean(Gender:SymCrossjk)", 
				"Gender:SymCrossjk:j:7;k:4 - mean(Gender:SymCrossjk)", "Gender:SymCrossjk:j:8;k:4 - mean(Gender:SymCrossjk)", 
				"Gender:SymCrossjk:j:6;k:5 - mean(Gender:SymCrossjk)", "Gender:SymCrossjk:j:7;k:5 - mean(Gender:SymCrossjk)", 
				"Gender:SymCrossjk:j:8;k:5 - mean(Gender:SymCrossjk)", "Gender:SymCrossjk:j:7;k:6 - mean(Gender:SymCrossjk)", 
				"Gender:SymCrossjk:j:8;k:6 - mean(Gender:SymCrossjk)", "Gender:SymCrossjk:j:8;k:7 - mean(Gender:SymCrossjk)", 
				"Gender:ASymCrossjkDkj:j:2;k:1 - mean(Gender:ASymCrossjkDkj)", "Gender:ASymCrossjkDkj:j:3;k:1 - mean(Gender:ASymCrossjkDkj)", 
				"Gender:ASymCrossjkDkj:j:4;k:1 - mean(Gender:ASymCrossjkDkj)", "Gender:ASymCrossjkDkj:j:5;k:1 - mean(Gender:ASymCrossjkDkj)", 
				"Gender:ASymCrossjkDkj:j:6;k:1 - mean(Gender:ASymCrossjkDkj)", "Gender:ASymCrossjkDkj:j:7;k:1 - mean(Gender:ASymCrossjkDkj)", 
				"Gender:ASymCrossjkDkj:j:8;k:1 - mean(Gender:ASymCrossjkDkj)", "Gender:ASymCrossjkDkj:j:3;k:2 - mean(Gender:ASymCrossjkDkj)", 
				"Gender:ASymCrossjkDkj:j:4;k:2 - mean(Gender:ASymCrossjkDkj)", "Gender:ASymCrossjkDkj:j:5;k:2 - mean(Gender:ASymCrossjkDkj)", 
				"Gender:ASymCrossjkDkj:j:6;k:2 - mean(Gender:ASymCrossjkDkj)", "Gender:ASymCrossjkDkj:j:7;k:2 - mean(Gender:ASymCrossjkDkj)", 
				"Gender:ASymCrossjkDkj:j:8;k:2 - mean(Gender:ASymCrossjkDkj)", "Gender:ASymCrossjkDkj:j:4;k:3 - mean(Gender:ASymCrossjkDkj)", 
				"Gender:ASymCrossjkDkj:j:5;k:3 - mean(Gender:ASymCrossjkDkj)", "Gender:ASymCrossjkDkj:j:6;k:3 - mean(Gender:ASymCrossjkDkj)", 
				"Gender:ASymCrossjkDkj:j:7;k:3 - mean(Gender:ASymCrossjkDkj)", "Gender:ASymCrossjkDkj:j:8;k:3 - mean(Gender:ASymCrossjkDkj)", 
				"Gender:ASymCrossjkDkj:j:5;k:4 - mean(Gender:ASymCrossjkDkj)", "Gender:ASymCrossjkDkj:j:6;k:4 - mean(Gender:ASymCrossjkDkj)", 
				"Gender:ASymCrossjkDkj:j:7;k:4 - mean(Gender:ASymCrossjkDkj)", "Gender:ASymCrossjkDkj:j:8;k:4 - mean(Gender:ASymCrossjkDkj)", 
				"Gender:ASymCrossjkDkj:j:6;k:5 - mean(Gender:ASymCrossjkDkj)", "Gender:ASymCrossjkDkj:j:7;k:5 - mean(Gender:ASymCrossjkDkj)", 
				"Gender:ASymCrossjkDkj:j:8;k:5 - mean(Gender:ASymCrossjkDkj)", "Gender:ASymCrossjkDkj:j:7;k:6 - mean(Gender:ASymCrossjkDkj)", 
				"Gender:ASymCrossjkDkj:j:8;k:6 - mean(Gender:ASymCrossjkDkj)", "Gender:ASymCrossjkDkj:j:8;k:7 - mean(Gender:ASymCrossjkDkj)")

new.labels <- c("v:B6;AJ", "v:129;AJ", "v:NOD;AJ", "v:NZO;AJ", "v:CAST;AJ", "v:PWK;AJ", "v:WSB;AJ", "v:129;B6", 
				"v:NOD;B6", "v:NZO;B6", "v:CAST;B6", "v:PWK;B6", "v:WSB;B6", "v:NOD;129", "v:NZO;129", "v:CAST;129", 
				"v:PWK;129", "v:WSB;129", "v:NZO;NOD", "v:CAST;NOD", "v:PWK;NOD", "v:WSB;NOD", "v:CAST;NZO", 
				"v:PWK;NZO", "v:WSB;NZO", "v:PWK;CAST", "v:WSB;CAST", "v:WSB;PWK", "w:B6;AJ", "w:129;AJ", "w:NOD;AJ", 
				"w:NZO;AJ", "w:CAST;AJ", "w:PWK;AJ", "w:WSB;AJ", "w:129;B6", "w:NOD;B6", "w:NZO;B6", "w:CAST;B6", "w:PWK;B6", 
				"w:WSB;B6", "w:NOD;129", "w:NZO;129", "w:CAST;129", "w:PWK;129", "w:WSB;129", "w:NZO;NOD", "w:CAST;NOD", 
				"w:PWK;NOD", "w:WSB;NOD", "w:CAST;NZO", "w:PWK;NZO", "w:WSB;NZO", "w:PWK;CAST", "w:WSB;CAST", "w:WSB;PWK")

plot.hpd(FluDi$AllDiallelObs[[1]]$centered.chains[,var.labels], names=new.labels, xlim=xlim, main="Sex/strainpair-specific") ## female-specific symmetric, asymmetric
## plot.hpd(FluDi$AllDiallelObs[[1]]$centered.chains[,c(134:161, 190:217)]) ## female-specific symmetric, asymmetric
abline(v=0, col="gray")

##par(mfrow=c(1,2))

##plot.hpd(FluDi$AllDiallelObs[[1]]$centered.chains[,218:229], cex.label=0.7) ## tau and sigma
##abline(v=0, col=gray)

##plot.hpd(FluDi$AllDiallelObs[[1]]$centered.chains[,6:57], cex.label=0.5) ## Batch
##abline(v=0, col="gray")
dev.off()

pdf(file.path(fname,"Heatmap_plots.pdf"), width=11, height=8.5)
par(mfrow=c(1,1))
## Plot Females, observed versus fitted
FluDi$AllDiallelObs[[1]]$TwoDiallelPlot(FemalePlot = TRUE,
MaleAgainstFemale = FALSE);
## Plot Males, observed versus fitted
FluDi$AllDiallelObs[[1]]$TwoDiallelPlot(FemalePlot = FALSE,
MaleAgainstFemale = FALSE);
## Plot Male Predicted versus Female Predicted.
FluDi$AllDiallelObs[[1]]$TwoDiallelPlot(FemalePlot = TRUE,
MaleAgainstFemale = TRUE);
## Plot Male Observed vresus Female Observed.
FluDi$AllDiallelObs[[1]]$TwoDiallelPlot(PlotOnlyObserved = TRUE)
dev.off()

pdf(file.path(fname,"Straw_plots.pdf"), width=11, height=8.5)
FluDi$AllDiallelObs[[1]]$PlotStrawPlot()
dev.off()

png(file.path(fname, "chains%03d.png"), width = 850, height = 1100)
plot(FluDi$AllDiallelObs[[1]]$raw.chains)
dev.off()
