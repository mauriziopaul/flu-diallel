library("treatmentResponseDiallel")
library("data.table")
library("coda")
library("xtable")

dat.D1 <- fread("/Users/maurizio/git/flu-diallel-paper/treatment-diallel/ValdarMachine/imputation-1000-D1-Diplo6/output/stackedTreatDiMIMPs.csv")
dat.D2 <- fread("/Users/maurizio/git/flu-diallel-paper/treatment-diallel/ValdarMachine/imputation-1000-D2-Diplo6/output/stackedTreatDiMIMPs.csv")
dat.D3 <- fread("/Users/maurizio/git/flu-diallel-paper/treatment-diallel/ValdarMachine/imputation-1000-D3-Diplo6/output/stackedTreatDiMIMPs.csv")
dat.D4 <- fread("/Users/maurizio/git/flu-diallel-paper/treatment-diallel/ValdarMachine/imputation-1000-D4-Diplo6/output/stackedTreatDiMIMPs.csv")
mx <- c("Batch:1:Mx1Diplo6:1", "Batch:1:Mx1Diplo6:2", "Batch:1:Mx1Diplo6:3", "Batch:1:Mx1Diplo6:4", "Batch:1:Mx1Diplo6:5", "Batch:1:Mx1Diplo6:6")

lower <- -0.75
upper <- 1.25
low <- -1
up <- 2
lwd <- 4

dat.mx <- as.data.frame(dat.D1)[,mx]
names(dat.mx) <- c("cast.cast", "cast.dom", "cast.mus", "dom.dom", "dom.mus", "mus.mus")
dat.mx$dominance.cast <- (dat.mx[,"cast.cast"]-dat.mx[,"cast.dom"])/(dat.mx[,"cast.cast"]-dat.mx[,"dom.dom"])
dat.mx$dominance.mus <- (dat.mx[,"mus.mus"]-dat.mx[,"dom.mus"])/(dat.mx[,"mus.mus"]-dat.mx[,"dom.dom"])
dat.mx$dominance.mus.v.cast <- (dat.mx[,"mus.mus"]-dat.mx[,"cast.mus"])/(dat.mx[,"mus.mus"]-dat.mx[,"cast.cast"])
dat.mx$mus.min.cast <- dat.mx$dominance.mus - dat.mx$dominance.cast
dat.mx$cast.min.mus <- dat.mx$dominance.cast - dat.mx$dominance.mus
dat.mx.D1 <- dat.mx
dd.D1 <- dat.mx$mus.min.cast

nrow(dat.mx[dat.mx$dominance.mus>=0.5,]) / nrow(dat.mx)
nrow(dat.mx[dat.mx$dominance.cast>=0.5,]) / nrow(dat.mx)
nrow(dat.mx[dat.mx$dominance.mus.v.cast>=0.5,]) / nrow(dat.mx)
nrow(dat.mx[dat.mx$dominance.mus<=0,]) / nrow(dat.mx)
nrow(dat.mx[dat.mx$dominance.cast<=0,]) / nrow(dat.mx)
nrow(dat.mx[dat.mx$dominance.mus.v.cast<=0,]) / nrow(dat.mx)
nrow(dat.mx[dat.mx$mus.min.cast<0,]) / nrow(dat.mx)

pdf("D1.dominance.cast.mus.pdf", width=6, height=5)
plot(density(dat.mx$dominance.mus, from=low, to=up), xlim=c(lower, upper), main="", col="blue", lwd=lwd,
	xlab="dominance index", yaxt="n", ylab="")
lines(density(dat.mx$dominance.mus.v.cast, from=low, to=up), col="orange", lwd=lwd)
lines(density(dat.mx$dominance.cast, from=low, to=up), col="green", lwd=lwd)
abline(v=c(-0.5, 0, 0.5, 1), col="grey", lty=2)
dev.off()

pdf("D1.dominance.diff.pdf", width=6, height=5)
plot(density(dat.mx$cast.min.mus, from=-1.25, to=1.25), xlim=c(-1, 1), main="", col="black", lwd=lwd,
	xlab="dominance difference index", yaxt="n", ylab="")
abline(v=0, col="black", lty=2)
abline(v=c(-1, 1), col="grey", lty=2)
dev.off()

dat.mx <- as.data.frame(dat.D2)[,mx]
names(dat.mx) <- c("cast.cast", "cast.dom", "cast.mus", "dom.dom", "dom.mus", "mus.mus")
dat.mx$dominance.cast <- (dat.mx[,"cast.cast"]-dat.mx[,"cast.dom"])/(dat.mx[,"cast.cast"]-dat.mx[,"dom.dom"])
dat.mx$dominance.mus <- (dat.mx[,"mus.mus"]-dat.mx[,"dom.mus"])/(dat.mx[,"mus.mus"]-dat.mx[,"dom.dom"])
dat.mx$dominance.mus.v.cast <- (dat.mx[,"mus.mus"]-dat.mx[,"cast.mus"])/(dat.mx[,"mus.mus"]-dat.mx[,"cast.cast"])
dat.mx$mus.min.cast <- dat.mx$dominance.mus - dat.mx$dominance.cast
dat.mx$cast.min.mus <- dat.mx$dominance.cast - dat.mx$dominance.mus
dat.mx.D2 <- dat.mx
dd.D2 <- dat.mx$mus.min.cast

nrow(dat.mx[dat.mx$dominance.mus>=0.5,]) / nrow(dat.mx)
nrow(dat.mx[dat.mx$dominance.cast>=0.5,]) / nrow(dat.mx)
nrow(dat.mx[dat.mx$dominance.mus.v.cast>=0.5,]) / nrow(dat.mx)
nrow(dat.mx[dat.mx$dominance.mus<=0,]) / nrow(dat.mx)
nrow(dat.mx[dat.mx$dominance.cast<=0,]) / nrow(dat.mx)
nrow(dat.mx[dat.mx$dominance.mus.v.cast<=0,]) / nrow(dat.mx)
nrow(dat.mx[dat.mx$mus.min.cast<0,]) / nrow(dat.mx)

pdf("D2.dominance.cast.mus.pdf", width=6, height=5)
plot(density(dat.mx$dominance.cast, from=low, to=up), xlim=c(lower, upper), main="", col="green", lwd=lwd,
	xlab="dominance index", yaxt="n", ylab="")
lines(density(dat.mx$dominance.mus.v.cast, from=low, to=up), col="orange", lwd=lwd)
lines(density(dat.mx$dominance.mus, from=low, to=up), col="blue", lwd=lwd)
abline(v=c(-0.5, 0, 0.5, 1), col="grey", lty=2)
dev.off()

pdf("D2.dominance.diff.pdf", width=6, height=5)
plot(density(dat.mx$mus.min.cast, from=-1.25, to=1.25), xlim=c(-1, 1), main="", col="black", lwd=lwd,
	xlab="dominance difference index", yaxt="n", ylab="")
abline(v=c(-1, 1), col="grey", lty=2)
abline(v=0, col="black", lty=2)
dev.off()

dat.mx <- as.data.frame(dat.D3)[,mx]
names(dat.mx) <- c("cast.cast", "cast.dom", "cast.mus", "dom.dom", "dom.mus", "mus.mus")
dat.mx$dominance.cast <- (dat.mx[,"cast.cast"]-dat.mx[,"cast.dom"])/(dat.mx[,"cast.cast"]-dat.mx[,"dom.dom"])
dat.mx$dominance.mus <- (dat.mx[,"mus.mus"]-dat.mx[,"dom.mus"])/(dat.mx[,"mus.mus"]-dat.mx[,"dom.dom"])
dat.mx$dominance.mus.v.cast <- (dat.mx[,"mus.mus"]-dat.mx[,"cast.mus"])/(dat.mx[,"mus.mus"]-dat.mx[,"cast.cast"])
dat.mx$mus.min.cast <- dat.mx$dominance.mus - dat.mx$dominance.cast
dat.mx$cast.min.mus <- dat.mx$dominance.cast - dat.mx$dominance.mus
dat.mx.D3 <- dat.mx
dd.D3 <- dat.mx$mus.min.cast

nrow(dat.mx[dat.mx$dominance.mus>=0.5,]) / nrow(dat.mx)
nrow(dat.mx[dat.mx$dominance.cast>=0.5,]) / nrow(dat.mx)
nrow(dat.mx[dat.mx$dominance.mus.v.cast>=0.5,]) / nrow(dat.mx)
nrow(dat.mx[dat.mx$dominance.mus<=0,]) / nrow(dat.mx)
nrow(dat.mx[dat.mx$dominance.cast<=0,]) / nrow(dat.mx)
nrow(dat.mx[dat.mx$dominance.mus.v.cast<=0,]) / nrow(dat.mx)
nrow(dat.mx[dat.mx$mus.min.cast<0,]) / nrow(dat.mx)

pdf("D3.dominance.cast.mus.pdf", width=6, height=5)
plot(density(dat.mx$dominance.cast, from=low, to=up), xlim=c(lower, upper), main="", col="green", lwd=lwd,
	xlab="dominance index", yaxt="n", ylab="")
lines(density(dat.mx$dominance.mus.v.cast, from=low, to=up), col="orange", lwd=lwd)
lines(density(dat.mx$dominance.mus, from=low, to=up), col="blue", lwd=lwd)
abline(v=c(-0.5, 0, 0.5, 1), col="grey", lty=2)
dev.off()

pdf("D3.dominance.diff.pdf", width=6, height=5)
plot(density(dat.mx$mus.min.cast, from=-1.25, to=1.25), xlim=c(-1, 1), main="", col="black", lwd=lwd,
	xlab="dominance difference index", yaxt="n", ylab="")
abline(v=c(-1, 1), col="grey", lty=2)
abline(v=0, col="black", lty=2)
dev.off()

dat.mx <- as.data.frame(dat.D4)[,mx]
names(dat.mx) <- c("cast.cast", "cast.dom", "cast.mus", "dom.dom", "dom.mus", "mus.mus")
dat.mx$dominance.cast <- (dat.mx[,"cast.cast"]-dat.mx[,"cast.dom"])/(dat.mx[,"cast.cast"]-dat.mx[,"dom.dom"])
dat.mx$dominance.mus <- (dat.mx[,"mus.mus"]-dat.mx[,"dom.mus"])/(dat.mx[,"mus.mus"]-dat.mx[,"dom.dom"])
dat.mx$dominance.mus.v.cast <- (dat.mx[,"mus.mus"]-dat.mx[,"cast.mus"])/(dat.mx[,"mus.mus"]-dat.mx[,"cast.cast"])
dat.mx$mus.min.cast <- dat.mx$dominance.mus - dat.mx$dominance.cast
dat.mx$cast.min.mus <- dat.mx$dominance.cast - dat.mx$dominance.mus
dat.mx.D4 <- dat.mx
dd.D4 <- dat.mx$mus.min.cast

nrow(dat.mx[dat.mx$dominance.mus>=0.5,]) / nrow(dat.mx)
nrow(dat.mx[dat.mx$dominance.cast>=0.5,]) / nrow(dat.mx)
nrow(dat.mx[dat.mx$dominance.mus.v.cast>=0.5,]) / nrow(dat.mx)
nrow(dat.mx[dat.mx$dominance.mus<=0,]) / nrow(dat.mx)
nrow(dat.mx[dat.mx$dominance.cast<=0,]) / nrow(dat.mx)
nrow(dat.mx[dat.mx$dominance.mus.v.cast<=0,]) / nrow(dat.mx)
nrow(dat.mx[dat.mx$mus.min.cast<0,]) / nrow(dat.mx)

pdf("D4.dominance.cast.mus.pdf", width=6, height=5)
plot(density(dat.mx$dominance.cast, from=low, to=up), xlim=c(lower, upper), main="", col="green", lwd=lwd,
	xlab="dominance index", yaxt="n", ylab="")
lines(density(dat.mx$dominance.mus.v.cast, from=low, to=up), col="orange", lwd=lwd)
lines(density(dat.mx$dominance.mus, from=low, to=up), col="blue", lwd=lwd)
abline(v=c(-0.5, 0, 0.5, 1), col="grey", lty=2)
dev.off()

pdf("D4.dominance.diff.pdf", width=6, height=5)
plot(density(dat.mx$mus.min.cast, from=-1.25, to=1.25), xlim=c(-1, 1), main="", col="black", lwd=lwd,
	xlab="dominance difference index", yaxt="n", ylab="")
abline(v=0, col="black", lty=2)
abline(v=c(-1, 1), col="grey", lty=2)
dev.off()

pdf("all.days.dominance.diff.pdf", width=6, height=5)
plot(density(dd.D4, from=-1.5, to=1.5), xlim=c(-1.25, 1.25), main="", col=rgb(0,0,0,maxColorValue=1), lwd=lwd,
	xlab="", yaxt="n", ylab="", cex.axis=1.5)
lines(density(dd.D3, from=-1.5, to=1.5), col=rgb(0,0,0,maxColorValue=1, alpha=0.5), lwd=lwd,
	xlab="dominance difference index", yaxt="n", ylab="")
lines(density(dd.D2, from=-1.5, to=1.5), col=rgb(0,0,0,maxColorValue=1, alpha=0.3), lwd=lwd,
	xlab="dominance difference index", yaxt="n", ylab="")
lines(density(dd.D1, from=-1.5, to=1.5), col=rgb(0,0,0,maxColorValue=1, alpha=0.1), lwd=lwd,
	xlab="dominance difference index", yaxt="n", ylab="")
abline(v=0, col="black", lty=1)
abline(v=c(-1, 1), col="grey", lty=2)
dev.off()

nrow(dat.mx[dat.mx$dominance.mus>=0.5,]) / nrow(dat.mx)
nrow(dat.mx[dat.mx$dominance.cast>=0.5,]) / nrow(dat.mx)
nrow(dat.mx[dat.mx$dominance.mus.v.cast>=0.5,]) / nrow(dat.mx)
nrow(dat.mx[dat.mx$dominance.mus<=0,]) / nrow(dat.mx)
nrow(dat.mx[dat.mx$dominance.cast<=0,]) / nrow(dat.mx)
nrow(dat.mx[dat.mx$dominance.mus.v.cast<=0,]) / nrow(dat.mx)
nrow(dat.mx[dat.mx$mus.min.cast<0,]) / nrow(dat.mx)

## COMPARISON

dat.mx <- as.data.frame(dat.D1)[,mx]
names(dat.mx) <- c("cast.cast", "cast.dom", "cast.mus", "dom.dom", "dom.mus", "mus.mus")
dat.mx$dominance.cast <- (dat.mx[,"cast.cast"]-dat.mx[,"cast.dom"])/(dat.mx[,"cast.cast"]-dat.mx[,"dom.dom"])
dat.mx$dominance.mus <- (dat.mx[,"mus.mus"]-dat.mx[,"dom.mus"])/(dat.mx[,"mus.mus"]-dat.mx[,"dom.dom"])
dat.mx$dominance.mus.v.cast <- (dat.mx[,"mus.mus"]-dat.mx[,"cast.mus"])/(dat.mx[,"mus.mus"]-dat.mx[,"cast.cast"])
dat.mx$mus.min.cast <- dat.mx$dominance.mus - dat.mx$dominance.cast
dat.mx$cast.min.mus <- dat.mx$dominance.cast - dat.mx$dominance.mus
dat.mx.D1 <- dat.mx
dd.D1 <- dat.mx$mus.min.cast

pdf("D1.dominance.cast.mus-comparison.pdf", width=6, height=10)
par(mfrow=c(2,1))
plot(density(dat.mx$dominance.mus, from=low, to=up), xlim=c(lower, upper), main="", col="blue", lwd=lwd,
     xlab="dominance index", yaxt="n", ylab="")
lines(density(dat.mx$dominance.mus.v.cast, from=low, to=up), col="orange", lwd=lwd)
lines(density(dat.mx$dominance.cast, from=low, to=up), col="green", lwd=lwd)
abline(v=c(-0.5, 0, 0.5, 1), col="grey", lty=2)
plot(density(dat.mx$dominance.mus), xlim=c(lower, upper), main="", col="blue", lwd=lwd,
     xlab="dominance index", yaxt="n", ylab="")
lines(density(dat.mx$dominance.mus.v.cast), col="orange", lwd=lwd)
lines(density(dat.mx$dominance.cast), col="green", lwd=lwd)
abline(v=c(-0.5, 0, 0.5, 1), col="grey", lty=2)
dev.off()

dat.mx <- as.data.frame(dat.D2)[,mx]
names(dat.mx) <- c("cast.cast", "cast.dom", "cast.mus", "dom.dom", "dom.mus", "mus.mus")
dat.mx$dominance.cast <- (dat.mx[,"cast.cast"]-dat.mx[,"cast.dom"])/(dat.mx[,"cast.cast"]-dat.mx[,"dom.dom"])
dat.mx$dominance.mus <- (dat.mx[,"mus.mus"]-dat.mx[,"dom.mus"])/(dat.mx[,"mus.mus"]-dat.mx[,"dom.dom"])
dat.mx$dominance.mus.v.cast <- (dat.mx[,"mus.mus"]-dat.mx[,"cast.mus"])/(dat.mx[,"mus.mus"]-dat.mx[,"cast.cast"])
dat.mx$mus.min.cast <- dat.mx$dominance.mus - dat.mx$dominance.cast
dat.mx$cast.min.mus <- dat.mx$dominance.cast - dat.mx$dominance.mus
dat.mx.D2 <- dat.mx
dd.D2 <- dat.mx$mus.min.cast

pdf("D2.dominance.cast.mus-comparison.pdf", width=6, height=10)
par(mfrow=c(2,1))
plot(density(dat.mx$dominance.mus, from=low, to=up), xlim=c(lower, upper), main="", col="blue", lwd=lwd,
     xlab="dominance index", yaxt="n", ylab="")
lines(density(dat.mx$dominance.mus.v.cast, from=low, to=up), col="orange", lwd=lwd)
lines(density(dat.mx$dominance.cast, from=low, to=up), col="green", lwd=lwd)
abline(v=c(-0.5, 0, 0.5, 1), col="grey", lty=2)
plot(density(dat.mx$dominance.mus), xlim=c(lower, upper), main="", col="blue", lwd=lwd,
     xlab="dominance index", yaxt="n", ylab="")
lines(density(dat.mx$dominance.mus.v.cast), col="orange", lwd=lwd)
lines(density(dat.mx$dominance.cast), col="green", lwd=lwd)
abline(v=c(-0.5, 0, 0.5, 1), col="grey", lty=2)
dev.off()

dat.mx <- as.data.frame(dat.D3)[,mx]
names(dat.mx) <- c("cast.cast", "cast.dom", "cast.mus", "dom.dom", "dom.mus", "mus.mus")
dat.mx$dominance.cast <- (dat.mx[,"cast.cast"]-dat.mx[,"cast.dom"])/(dat.mx[,"cast.cast"]-dat.mx[,"dom.dom"])
dat.mx$dominance.mus <- (dat.mx[,"mus.mus"]-dat.mx[,"dom.mus"])/(dat.mx[,"mus.mus"]-dat.mx[,"dom.dom"])
dat.mx$dominance.mus.v.cast <- (dat.mx[,"mus.mus"]-dat.mx[,"cast.mus"])/(dat.mx[,"mus.mus"]-dat.mx[,"cast.cast"])
dat.mx$mus.min.cast <- dat.mx$dominance.mus - dat.mx$dominance.cast
dat.mx$cast.min.mus <- dat.mx$dominance.cast - dat.mx$dominance.mus
dat.mx.D3 <- dat.mx
dd.D3 <- dat.mx$mus.min.cast

pdf("D3.dominance.cast.mus-comparison.pdf", width=6, height=10)
par(mfrow=c(2,1))
plot(density(dat.mx$dominance.mus, from=low, to=up), xlim=c(lower, upper), main="", col="blue", lwd=lwd,
     xlab="dominance index", yaxt="n", ylab="")
lines(density(dat.mx$dominance.mus.v.cast, from=low, to=up), col="orange", lwd=lwd)
lines(density(dat.mx$dominance.cast, from=low, to=up), col="green", lwd=lwd)
abline(v=c(-0.5, 0, 0.5, 1), col="grey", lty=2)
plot(density(dat.mx$dominance.mus), xlim=c(lower, upper), main="", col="blue", lwd=lwd,
     xlab="dominance index", yaxt="n", ylab="")
lines(density(dat.mx$dominance.mus.v.cast), col="orange", lwd=lwd)
lines(density(dat.mx$dominance.cast), col="green", lwd=lwd)
abline(v=c(-0.5, 0, 0.5, 1), col="grey", lty=2)
dev.off()

dat.mx <- as.data.frame(dat.D4)[,mx]
names(dat.mx) <- c("cast.cast", "cast.dom", "cast.mus", "dom.dom", "dom.mus", "mus.mus")
dat.mx$dominance.cast <- (dat.mx[,"cast.cast"]-dat.mx[,"cast.dom"])/(dat.mx[,"cast.cast"]-dat.mx[,"dom.dom"])
dat.mx$dominance.mus <- (dat.mx[,"mus.mus"]-dat.mx[,"dom.mus"])/(dat.mx[,"mus.mus"]-dat.mx[,"dom.dom"])
dat.mx$dominance.mus.v.cast <- (dat.mx[,"mus.mus"]-dat.mx[,"cast.mus"])/(dat.mx[,"mus.mus"]-dat.mx[,"cast.cast"])
dat.mx$mus.min.cast <- dat.mx$dominance.mus - dat.mx$dominance.cast
dat.mx$cast.min.mus <- dat.mx$dominance.cast - dat.mx$dominance.mus
dat.mx.D4 <- dat.mx
dd.D4 <- dat.mx$mus.min.cast

pdf("D4.dominance.cast.mus-comparison.pdf", width=6, height=10)
par(mfrow=c(2,1))
plot(density(dat.mx$dominance.mus, from=low, to=up), xlim=c(lower, upper), main="", col="blue", lwd=lwd,
     xlab="dominance index", yaxt="n", ylab="")
lines(density(dat.mx$dominance.mus.v.cast, from=low, to=up), col="orange", lwd=lwd)
lines(density(dat.mx$dominance.cast, from=low, to=up), col="green", lwd=lwd)
abline(v=c(-0.5, 0, 0.5, 1), col="grey", lty=2)
plot(density(dat.mx$dominance.mus), xlim=c(lower, upper), main="", col="blue", lwd=lwd,
     xlab="dominance index", yaxt="n", ylab="")
lines(density(dat.mx$dominance.mus.v.cast), col="orange", lwd=lwd)
lines(density(dat.mx$dominance.cast), col="green", lwd=lwd)
abline(v=c(-0.5, 0, 0.5, 1), col="grey", lty=2)
dev.off()

## Composite Image
pdf("D3.D4.dominance.cast.mus.pdf", width=6, height=12)
par(mfrow=c(3,1))
dat.mx <- as.data.frame(dat.D3)[,mx]
names(dat.mx) <- c("cast.cast", "cast.dom", "cast.mus", "dom.dom", "dom.mus", "mus.mus")
dat.mx$dominance.cast <- (dat.mx[,"cast.cast"]-dat.mx[,"cast.dom"])/(dat.mx[,"cast.cast"]-dat.mx[,"dom.dom"])
dat.mx$dominance.mus <- (dat.mx[,"mus.mus"]-dat.mx[,"dom.mus"])/(dat.mx[,"mus.mus"]-dat.mx[,"dom.dom"])
dat.mx$dominance.mus.v.cast <- (dat.mx[,"mus.mus"]-dat.mx[,"cast.mus"])/(dat.mx[,"mus.mus"]-dat.mx[,"cast.cast"])
dat.mx$mus.min.cast <- dat.mx$dominance.mus - dat.mx$dominance.cast
dat.mx$cast.min.mus <- dat.mx$dominance.cast - dat.mx$dominance.mus
dd.D3 <- dat.mx$mus.min.cast

par(mar=c(0,1,5,1))
plot(density(dat.mx$dominance.cast, from=low, to=up), xlim=c(lower, upper), main="", col="green", lwd=lwd,
     yaxt="n", ylab="", xaxt="n", xlab="", cex.axis=2)
lines(density(dat.mx$dominance.mus.v.cast, from=low, to=up), col="orange", lwd=lwd)
lines(density(dat.mx$dominance.mus, from=low, to=up), col="blue", lwd=lwd)
abline(v=c(-0.5, 0, 0.5, 1), col="grey", lty=2)

dat.mx <- as.data.frame(dat.D4)[,mx]
names(dat.mx) <- c("cast.cast", "cast.dom", "cast.mus", "dom.dom", "dom.mus", "mus.mus")
dat.mx$dominance.cast <- (dat.mx[,"cast.cast"]-dat.mx[,"cast.dom"])/(dat.mx[,"cast.cast"]-dat.mx[,"dom.dom"])
dat.mx$dominance.mus <- (dat.mx[,"mus.mus"]-dat.mx[,"dom.mus"])/(dat.mx[,"mus.mus"]-dat.mx[,"dom.dom"])
dat.mx$dominance.mus.v.cast <- (dat.mx[,"mus.mus"]-dat.mx[,"cast.mus"])/(dat.mx[,"mus.mus"]-dat.mx[,"cast.cast"])
dat.mx$mus.min.cast <- dat.mx$dominance.mus - dat.mx$dominance.cast
dat.mx$cast.min.mus <- dat.mx$dominance.cast - dat.mx$dominance.mus
dat.mx.d4 <- dat.mx
dd.D4 <- dat.mx$mus.min.cast

par(mar=c(5,1,0,1))
plot(density(dat.mx$dominance.cast, from=low, to=up), xlim=c(lower, upper), main="", col="green", lwd=lwd,
     xlab="", yaxt="n", ylab="", cex.axis=2)
lines(density(dat.mx$dominance.mus.v.cast, from=low, to=up), col="orange", lwd=lwd)
lines(density(dat.mx$dominance.mus, from=low, to=up), col="blue", lwd=lwd)
abline(v=c(-0.5, 0, 0.5, 1), col="grey", lty=2)

plot(density(dd.D4, from=-1, to=2), xlim=c(-1.25, 1.25), main="", col=rgb(0,0,0,maxColorValue=1), lwd=lwd,
     xlab="", yaxt="n", ylab="", cex.axis=2)
lines(density(dd.D3, from=-1.5, to=1.5), col=rgb(0,0,0,maxColorValue=1, alpha=0.5), lwd=lwd,
      xlab="dominance difference index", yaxt="n", ylab="")
lines(density(dd.D2, from=-1.5, to=1.5), col=rgb(0,0,0,maxColorValue=1, alpha=0.3), lwd=lwd,
      xlab="dominance difference index", yaxt="n", ylab="")
lines(density(dd.D1, from=-1.5, to=1.5), col=rgb(0,0,0,maxColorValue=1, alpha=0.1), lwd=lwd,
      xlab="dominance difference index", yaxt="n", ylab="")
abline(v=0, col="black", lty=1)
abline(v=c(-1, 1), col="grey", lty=2)
dev.off()

density(dd.D1, from=-1.5, to=1.5)
d1 <- density(dd.D1,from=-1.5, to=1.5)
d2 <- density(dd.D2,from=-1.5, to=1.5)
d3 <- density(dd.D3,from=-1.5, to=1.5)
d4 <- density(dd.D4,from=-1.5, to=1.5)
plot(d1$y~d1$x, type="l", ylim=c(0,1))
lines(d2$y~d2$x)
lines(d3$y~d3$x)
lines(d4$y~d4$x)

# posterior mode, dominance difference
d1$x[d1$y==max(d1$y)]
d2$x[d2$y==max(d2$y)]
d3$x[d3$y==max(d3$y)]
d4$x[d4$y==max(d4$y)]

HPDinterval(as.mcmc(dd.D1), prob=0.5)[,1:2]; median(dd.D1)
HPDinterval(as.mcmc(dd.D2), prob=0.5)[,1:2]; median(dd.D2)
HPDinterval(as.mcmc(dd.D3), prob=0.5)[,1:2]; median(dd.D3)
HPDinterval(as.mcmc(dd.D4), prob=0.5)[,1:2]; median(dd.D4)

HPDinterval(as.mcmc(dat.mx.D1$dominance.cast), prob=0.5)[,1:2]
HPDinterval(as.mcmc(dat.mx.D2$dominance.cast), prob=0.5)[,1:2]
HPDinterval(as.mcmc(dat.mx.D3$dominance.cast), prob=0.5)[,1:2]
HPDinterval(as.mcmc(dat.mx.D4$dominance.cast), prob=0.5)[,1:2]
HPDinterval(as.mcmc(dat.mx.D1$dominance.cast), prob=0.8)[,1:2]
HPDinterval(as.mcmc(dat.mx.D2$dominance.cast), prob=0.8)[,1:2]
HPDinterval(as.mcmc(dat.mx.D3$dominance.cast), prob=0.8)[,1:2]
HPDinterval(as.mcmc(dat.mx.D4$dominance.cast), prob=0.8)[,1:2]
median(dat.mx.D1$dominance.cast)
median(dat.mx.D2$dominance.cast)
median(dat.mx.D3$dominance.cast)
median(dat.mx.D4$dominance.cast)

d1.cast <- density(dat.mx.D1$dominance.cast,from=-1, to=2)
d2.cast <- density(dat.mx.D2$dominance.cast,from=-1, to=2)
d3.cast <- density(dat.mx.D3$dominance.cast,from=-1, to=2)
d4.cast <- density(dat.mx.D4$dominance.cast,from=-1, to=2)
plot(d1.cast$y~d1.cast$x, type="l", ylim=c(0,2))
lines(d2.cast$y~d2.cast$x)
lines(d3.cast$y~d3.cast$x)
lines(d4.cast$y~d4.cast$x)

## Posterior mode
d1.cast$x[d1.cast$y==max(d1.cast$y)]
d2.cast$x[d2.cast$y==max(d2.cast$y)]
d3.cast$x[d3.cast$y==max(d3.cast$y)]
d4.cast$x[d4.cast$y==max(d4.cast$y)]

HPDinterval(as.mcmc(dat.mx.D1$dominance.mus), prob=0.5)[,1:2]
HPDinterval(as.mcmc(dat.mx.D2$dominance.mus), prob=0.5)[,1:2]
HPDinterval(as.mcmc(dat.mx.D3$dominance.mus), prob=0.5)[,1:2]
HPDinterval(as.mcmc(dat.mx.D4$dominance.mus), prob=0.5)[,1:2]
HPDinterval(as.mcmc(dat.mx.D1$dominance.mus), prob=0.8)[,1:2]
HPDinterval(as.mcmc(dat.mx.D2$dominance.mus), prob=0.8)[,1:2]
HPDinterval(as.mcmc(dat.mx.D3$dominance.mus), prob=0.8)[,1:2]
HPDinterval(as.mcmc(dat.mx.D4$dominance.mus), prob=0.8)[,1:2]
median(dat.mx.D1$dominance.mus)
median(dat.mx.D2$dominance.mus)
median(dat.mx.D3$dominance.mus)
median(dat.mx.D4$dominance.mus)

d1.mus <- density(dat.mx.D1$dominance.mus,from=-1, to=2)
d2.mus <- density(dat.mx.D2$dominance.mus,from=-1, to=2)
d3.mus <- density(dat.mx.D3$dominance.mus,from=-1, to=2)
d4.mus <- density(dat.mx.D4$dominance.mus,from=-1, to=2)
plot(d1.mus$y~d1.mus$x, type="l", ylim=c(0,1.5))
lines(d2.mus$y~d2.mus$x)
lines(d3.mus$y~d3.mus$x)
lines(d4.mus$y~d4.mus$x)

## Posterior mode
d1.mus$x[d1.mus$y==max(d1.mus$y)]
d2.mus$x[d2.mus$y==max(d2.mus$y)]
d3.mus$x[d3.mus$y==max(d3.mus$y)]
d4.mus$x[d4.mus$y==max(d4.mus$y)]


HPDinterval(as.mcmc(dat.mx.D1$dominance.mus.v.cast), prob=0.5)[,1:2]
HPDinterval(as.mcmc(dat.mx.D2$dominance.mus.v.cast), prob=0.5)[,1:2]
HPDinterval(as.mcmc(dat.mx.D3$dominance.mus.v.cast), prob=0.5)[,1:2]
HPDinterval(as.mcmc(dat.mx.D4$dominance.mus.v.cast), prob=0.5)[,1:2]
HPDinterval(as.mcmc(dat.mx.D1$dominance.mus.v.cast), prob=0.8)[,1:2]
HPDinterval(as.mcmc(dat.mx.D2$dominance.mus.v.cast), prob=0.8)[,1:2]
HPDinterval(as.mcmc(dat.mx.D3$dominance.mus.v.cast), prob=0.8)[,1:2]
HPDinterval(as.mcmc(dat.mx.D4$dominance.mus.v.cast), prob=0.8)[,1:2]
median(dat.mx.D1$dominance.mus.v.cast)
median(dat.mx.D2$dominance.mus.v.cast)
median(dat.mx.D3$dominance.mus.v.cast)
median(dat.mx.D4$dominance.mus.v.cast)

d1.mus.v.cast <- density(dat.mx.D1$dominance.mus.v.cast,from=-1, to=2)
d2.mus.v.cast <- density(dat.mx.D2$dominance.mus.v.cast,from=-1, to=2)
d3.mus.v.cast <- density(dat.mx.D3$dominance.mus.v.cast,from=-1, to=2)
d4.mus.v.cast <- density(dat.mx.D4$dominance.mus.v.cast,from=-1, to=2)
plot(d1.mus.v.cast$y~d1.mus.v.cast$x, type="l", ylim=c(0,1.5))
lines(d2.mus.v.cast$y~d2.mus.v.cast$x)
lines(d3.mus.v.cast$y~d3.mus.v.cast$x)
lines(d4.mus.v.cast$y~d4.mus.v.cast$x)

# posterior mode
d1.mus.v.cast$x[d1.mus.v.cast$y==max(d1.mus.v.cast$y)]
d2.mus.v.cast$x[d2.mus.v.cast$y==max(d2.mus.v.cast$y)]
d3.mus.v.cast$x[d3.mus.v.cast$y==max(d3.mus.v.cast$y)]
d4.mus.v.cast$x[d4.mus.v.cast$y==max(d4.mus.v.cast$y)]

# csv to tex
dom.table <- read.csv("dominance-results.csv")
xtable(dom.table, digits=3)
