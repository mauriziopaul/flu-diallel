#!/usr/local/bin/Rscript

##---------------------------------------------------------------------------------------------------------------------
## Title: MakePostPlots Script
## Author: Paul L. Maurizio
## Institution: UNC-Chapel Hill
## Program: Bioinformatics and Computational Biology Ph.D. Curriculum
## Advisors: Will Valdar, Mark Heise
## Date Created: 2016-05-09
## Date Updated: 2017-01-05
##
## To Do:
##---------------------------------------------------------------------------------------------------------------------

library(cmdline)
library(configfile)
library(treatmentResponseDiallel)
library(tools)

confile <- cmdline.string("configfile")
savedirs <- cmdline.strings("savedirs")
config <- read.configfile(confile)

print(savedirs)

#dir.create(savedir, recursive=TRUE, showWarnings=FALSE)
source("loadConfig.R")

##-----------------------------------------------
## Analyze Stacked
##-----------------------------------------------

## PLOTS

# savedir <- "./All/pct_D1/Mx1diplo/starting_weight/MP/two-to-one/NumMatches_002/"
# savedir <- "./All/pct_D2/Mx1diplo/starting_weight/MP/two-to-one/NumMatches_002/"
# savedir <- "./All/pct_D3/Mx1diplo/starting_weight/MP/two-to-one/NumMatches_002/"
# savedir <- "./All/pct_D4/Mx1Diplo/starting_weight/MP/two-to-one/NumMatches_002/"

# savedirs <- c("./All/pct_D1/Mx1Diplo/starting_weight/MP/two-to-one/NumMatches_002/",
#  "./All/pct_D2/Mx1Diplo/starting_weight/MP/two-to-one/NumMatches_002/",
#  "./All/pct_D3/Mx1Diplo/starting_weight/MP/two-to-one/NumMatches_002/",
#  "./All/pct_D4/Mx1Diplo/starting_weight/MP/two-to-one/NumMatches_002/")

#savedirs <- c(
#		"./imputation-1000-D1-Diplo6",
#		"./imputation-1000-D2-Diplo6",
#		"./imputation-1000-D3-Diplo6",
#		"./imputation-1000-D4-Diplo6"
#		)

plot.mx.effects <- function(savedirs, batchIDs, batchNames, muID, batchColors="grey", ...){
	rowdim <- length(batchIDs)
	coldim <- length(savedirs)
	dat <- matrix(data=rep(NA, rowdim*(coldim+1)), nrow=rowdim, ncol=coldim+1)
	dat[,1] <- 100
	colnames(dat) <- c("D0", paste("D", c(1:coldim), sep=""))
	rownames(dat) <- batchNames
	for(i in c(1:coldim)){ # i index the days including D0
		stacked.matches <- as.mcmc(read.csv(file.path(savedirs[i], "output/stackedTreatDiMIMPs.csv"),
			check.names=FALSE))
		mu <- stacked.matches[,muID]
		mx <- stacked.matches[,batchIDs]
#		rowmeans <- rowMeans(as.matrix(mx))
#		mx <- mx - rowmeans
		mx <- mx + mu + 100
		colmeans <- colMeans(mx)
		dat[,i+1] <- colmeans
		colnames(mx) <- batchNames
		mx.unordered <- mx
		mx <- mx[,order(colmeans, decreasing=TRUE)]
		
		pdf(file.path(savedirs[i], paste0("output/mx1_timecourse_on_", colnames(dat)[i+1], ".pdf")), width=6, height=6)
		plot(dat[1,], ylim=c(85,101), type="b", col="white", lwd=2,
			ylab="pct starting wt", xlab="DPI", xlim=c(1,6.5), xaxt="n");
		for(j in c(1:rowdim)){ #j index the batches
			lines(dat[j,], type="b", col=batchColors[j], lwd=2);
			text(x=5.5, y=dat[i,], labels=names(dat[j,i+1]), pos=4, cex=0.9, col=batchColors[j])
		}
		# yheight <- c(97, 95, 94, 93, 91, 87)
		# text(x=i+1.5, y=dat[i,], labels=names(dat), pos=4, cex=0.9, col=batchColors)
		#text(x=5.5, y=86+(rank(dat[,i+1])), labels=names(dat[,i+1]), pos=4, cex=0.9, col=batchColors)
		#text(x=i+1.5, y=dat[i,], labels=names(dat), pos=4, cex=0.9, col=batchColors)
		axis(1, at=c(1:5), labels=c(0:coldim))
		abline(h=100, lty=2)
		dev.off()

		pdf(file.path(savedirs[i], "output", 
			paste("mx1_effects_on_", colnames(dat)[i+1],".pdf", sep="")), 
			width=4, height=6)
		plot.hpd(mx, ylab="", 
			xlab=paste("pct wt loss on ", as.character(colnames(dat)[i+1]), sep=""), 
			xlim=c(85,101))
		abline(v=100, col="grey")
		dev.off()

		pdf(file.path(savedirs[i], "output", 
			paste("mx1_unordered_effects_on_", colnames(dat)[i+1],".pdf", sep="")), 
			width=4, height=6)
		plot.hpd(mx.unordered, ylab="", 
			xlab=paste("pct wt loss on ", as.character(colnames(dat)[i+1]), sep=""), 
			xlim=c(85,101))
		abline(v=100, col="grey")
		dev.off()

	}

}

## Mx1 `Batch Effect` (+ mu)
## from stacked.matches

batchIDs <- c(
	"Batch:1:Mx1Diplo6:1",
	"Batch:1:Mx1Diplo6:2",
	"Batch:1:Mx1Diplo6:3",
	"Batch:1:Mx1Diplo6:4",
	"Batch:1:Mx1Diplo6:5",
	"Batch:1:Mx1Diplo6:6")

#[1] "Mx1Diplo6:dom_dom=4"   "Mx1Diplo6:dom_mus=5"   "Mx1Diplo6:cast_dom=2"
#[4] "Mx1Diplo6:mus_mus=6"   "Mx1Diplo6:cast_mus=3"  "Mx1Diplo6:cast_cast=1"
# [1] "Block:12=11" "Block:25=20" "Block:23=19" "Block:56=35" "Block:8=7"
# [6] "Block:15=13" "Block:64=38" "Block:14=12" "Block:43=32" "Block:84=49"
#[11] "Block:51=33" "Block:57=36" "Block:82=48" "Block:59=37" "Block:91=51"
#[16] "Block:69=40" "Block:90=50" "Block:32=25" "Block:5=4"   "Block:21=17"
#[21] "Block:1=1"   "Block:7=6"   "Block:19=15" "Block:33=26" "Block:34=27"
#[26] "Block:36=28" "Block:30=23" "Block:27=21" "Block:40=31" "Block:11=10"
#[31] "Block:20=16" "Block:2=2"   "Block:6=5"   "Block:38=30" "Block:31=24"
#[36] "Block:37=29" "Block:22=18" "Block:92=52" "Block:77=46" "Block:74=44"
#[41] "Block:18=14" "Block:72=43" "Block:75=45" "Block:3=3"   "Block:9=8"
#[46] "Block:29=22" "Block:71=42" "Block:80=47" "Block:53=34" "Block:66=39"
#[51] "Block:70=41" "Block:10=9"
# Blocks might differ

batchNames <- c(
	"cast. x cast.",
	"cast. x dom.",
	"cast. x mus.",
	"dom. x dom.",
	"dom. x mus.",
	"mus. x mus.")

batchColors <- c(
			rgb(0,160,0, maxColorValue=255),
			"black",
			"black",
			rgb(128,128,128, maxColorValue=255),
			"black",
			"blue"
#			rgb(0,0,240, maxColorValue=255)
#			rgb(240,0,0, maxColorValue=255)
		)

muID <- "Mu"
dat <- plot.mx.effects(savedirs, batchIDs, batchNames, muID, batchColors)

cat("Finished plotting! \n")

# stacked.posterior <- read.delim2(file.path(savedir[1], "output/stackedPostTreatDiMIMPs.csv"), 
# 	check.names=FALSE, sep=";", dec=",")

# add <- c(
# 	"additive:AJ",
# 	"additive:B6",
# 	"additive:129",
# 	"additive:NOD",
# 	"additive:NZO",
# 	"additive:CAST",
# 	"additive:PWK",
# 	"additive:WSB")

## epis

# v.class <- c(
# 	'v:B6;AJ','v:129;AJ','v:NOD;AJ','v:NZO;AJ','v:CAST;AJ',
# 	'v:PWK;AJ','v:WSB;AJ','v:129;B6','v:NOD;B6','v:NZO;B6',
# 	'v:CAST;B6','v:PWK;B6','v:WSB;B6','v:NOD;129','v:NZO;129',
# 	'v:CAST;129','v:PWK;129','v:WSB;129','v:NZO;NOD','v:CAST;NOD',
# 	'v:PWK;NOD','v:WSB;NOD','v:CAST;NZO','v:PWK;NZO','v:WSB;NZO',
# 	'v:PWK;CAST','v:WSB;CAST','v:WSB;PWK')
# v1 <- stacked.matches[,v.class]
# rowmeans <- rowMeans(as.matrix(v1))
# v2 <- v1 - rowmeans

# par(mfrow=c(1,2))
# xlim <- c(-10,10)
# plot.hpd(v1, xlim=xlim); abline(v=0, col="grey")
# plot.hpd(v2, xlim=xlim); abline(v=0, col="grey")

## Gridplot 
## from stacked.postpred

# stacked.postpred <- as.mcmc(read.csv(file.path(savedir[1], "output/stackedPostPredTreatDiMIMPs.csv"), 
# 	check.names=FALSE))

# pp <- stacked.postpred
# pp.S <- c(rep(0,64),rep(1,64))
# pp.i <- rep(seq(1,8,1),16)
# pp.k <- rep(
# 	c(rep(1,8), rep(2,8), rep(3,8), rep(4,8),
# 		rep(5,8), rep(6,8), rep(7,8), rep(8,8)),
# 	2)

# par(mfrow=c(8,8))
# for(S in c(0:1)){
# 	for(i in c(1:8)){
# 		for(k in c(1:8)){
# 			plot(1,1,main=paste(S, j, k))
# 		}
# 	}
# }

# par(mfrow=c(1,4))
# xlim <- c(-20, 5)

# plot.hpd(pp[,c(1:32)], xlim=xlim)
# plot.hpd(pp[,c(65:96)], xlim=xlim, add=TRUE, bg="grey")
# abline(v=0, col="grey")
# plot.hpd(pp[,c(33:64)], xlim=xlim)
# abline(v=0, col="grey")
# plot.hpd(pp[,c(65:96)], xlim=xlim)
# abline(v=0, col="grey")
# plot.hpd(pp[,c(97:128)], xlim=xlim)
# abline(v=0, col="grey")
