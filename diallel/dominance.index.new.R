library("treatmentResponseDiallel")
library("data.table")

dat.D1 <- fread("./imputation-1000-D1-Diplo6/output/stackedTreatDiMIMPs.csv")
dat.D2 <- fread("./imputation-1000-D2-Diplo6/output/stackedTreatDiMIMPs.csv")
dat.D3 <- fread("./imputation-1000-D3-Diplo6/output/stackedTreatDiMIMPs.csv")
dat.D4 <- fread("./imputation-1000-D4-Diplo6/output/stackedTreatDiMIMPs.csv")
mx <- c("Batch:1:Mx1Diplo6:1", "Batch:1:Mx1Diplo6:2", "Batch:1:Mx1Diplo6:3", "Batch:1:Mx1Diplo6:4", "Batch:1:Mx1Diplo6:5", "Batch:1:Mx1Diplo6:6")

lower <- -0.75
upper <- 1.25
low <- -1
up <- 1.5
lwd <- 3

dat.mx <- as.data.frame(dat.D1)[,mx]
names(dat.mx) <- c("cast.cast", "cast.dom", "cast.mus", "dom.dom", "dom.mus", "mus.mus")
dat.mx$dominance.cast <- (dat.mx[,"cast.cast"]-dat.mx[,"cast.dom"])/(dat.mx[,"cast.cast"]-dat.mx[,"dom.dom"])
dat.mx$dominance.mus <- (dat.mx[,"mus.mus"]-dat.mx[,"dom.mus"])/(dat.mx[,"mus.mus"]-dat.mx[,"dom.dom"])
dat.mx$dominance.mus.v.cast <- (dat.mx[,"mus.mus"]-dat.mx[,"cast.mus"])/(dat.mx[,"mus.mus"]-dat.mx[,"cast.cast"])
dat.mx$mus.min.cast <- dat.mx$dominance.mus - dat.mx$dominance.cast
dat.mx$cast.min.mus <- dat.mx$dominance.cast - dat.mx$dominance.mus
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

lwd <- 4
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

## Composite Image

pdf("D3.D4.dominance.cast.mus.pdf", width=6, height=14)
par(mfrow=c(3,1))

dat.mx <- as.data.frame(dat.D3)[,mx]
names(dat.mx) <- c("cast.cast", "cast.dom", "cast.mus", "dom.dom", "dom.mus", "mus.mus")
dat.mx$dominance.cast <- (dat.mx[,"cast.cast"]-dat.mx[,"cast.dom"])/(dat.mx[,"cast.cast"]-dat.mx[,"dom.dom"])
dat.mx$dominance.mus <- (dat.mx[,"mus.mus"]-dat.mx[,"dom.mus"])/(dat.mx[,"mus.mus"]-dat.mx[,"dom.dom"])
dat.mx$dominance.mus.v.cast <- (dat.mx[,"mus.mus"]-dat.mx[,"cast.mus"])/(dat.mx[,"mus.mus"]-dat.mx[,"cast.cast"])
dat.mx$mus.min.cast <- dat.mx$dominance.mus - dat.mx$dominance.cast
dat.mx$cast.min.mus <- dat.mx$dominance.cast - dat.mx$dominance.mus
dd.D3 <- dat.mx$mus.min.cast

par(mar=c(0,1,7,1))
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
dd.D4 <- dat.mx$mus.min.cast

par(mar=c(7,1,0,1))
plot(density(dat.mx$dominance.cast, from=low, to=up), xlim=c(lower, upper), main="", col="green", lwd=lwd,
     xlab="", yaxt="n", ylab="", cex.axis=2)
lines(density(dat.mx$dominance.mus.v.cast, from=low, to=up), col="orange", lwd=lwd)
lines(density(dat.mx$dominance.mus, from=low, to=up), col="blue", lwd=lwd)
abline(v=c(-0.5, 0, 0.5, 1), col="grey", lty=2)

plot(density(dd.D4, from=-1.5, to=1.5), xlim=c(-1.25, 1.25), main="", col=rgb(0,0,0,maxColorValue=1), lwd=lwd,
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
