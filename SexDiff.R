#!/usr/bin/Rscript

data <- read.csv("Flu_Diallel_v007.csv")
pct_D4 <- 100*data$D4/data$D0
data <- cbind(data, pct_D4)
model <- aggregate(pct_D4 ~ Strain + Sex + Trt, data=data, FUN="mean")
## 246 rows becomes 122 rows. We lost one pair, FD or CASTxNOD

males <- subset(model, Sex=="M")
females <- subset(model, Sex=="F")

## Reshape / merge so that male and female weights are in adjacent columns
m1 <- merge(males, females, by.x = c("Strain","Trt"), by.y = c("Strain","Trt"))
m1 <- subset(m1, Trt=="FLU")
pct_m_minus_f <- pct_m_minus_f <- m1$pct_D4.x - m1$pct_D4.y
m1 <- cbind(m1, pct_m_minus_f)
m1.sorted <- m1[order(pct_m_minus_f),]

pdf("Sex_Diff_Distribution.pdf")
plot(pct_m_minus_f ~ c(1:length(m1.sorted$Strain)), data=m1.sorted, 
	xaxt='n', pch=16, cex=0.7, ylab="D4 weight loss, M-F (% D0 weight)", xlab="")
abline(h=c(0, -10, 10), lty=2, col="gray")
dev.off()

pdf("Sex_Diff_Distribution_Hist.pdf")
hist(m1.sorted$pct_m_minus_f, xlab="D4 weight loss, M - F (% D0 weight)", main=NULL)
plot(density(m1.sorted$pct_m_minus_f), xlab="D4 weight loss, M - F (% D0 weight)", main="", yaxt='n')
abline(v=c(0,10, -10, -20, 20), col="gray", lty=2)
dev.off()

## List the crosses where the difference is more than -10 or +10
m1.tails <- subset(m1, abs(pct_m_minus_f) > 10)
strain.extremes <- c("DA", "DB", "DC", "EB", "GD", "GF", "HB", "HC")

## Plot extremes

pdf("sex_diff_extremes.pdf")
par(mfrow=c(2,2))
strain.list <- strain.extremes
for(i in strain.list){
	dat <- subset(data, Strain==i)
	dat <- subset(dat, (!is.na(pct_D4)))
	dat$Trt <- factor(dat$Trt)
	dat$Sex <- factor(dat$Sex)
	dat$Strain <- factor(dat$Strain)
	dat$Block <- factor(dat$Block)
	boxplot(pct_D4 ~ Strain*Trt*Sex + Block, data=dat, cex.axis=0.7, border="gray", boxwex=0.2, 
			ylim=c(70, 105), main=i)
	beeswarm(pct_D4 ~ Strain*Trt*Sex + Block, data=dat, add=TRUE, pch=16)
	dat <- subset(dat, Trt=="FLU")
	boxplot(pct_D4 ~ Strain*Sex + Block, data=dat, cex.axis=0.7, border="gray", boxwex=0.2, 
			ylim=c(70,105), main=i)
	beeswarm(pct_D4 ~ Strain*Sex + Block, data=dat, add=TRUE, pch=16)
}
dev.off()

plot(pct_D4.x ~ pct_D4.y, data=m1, ylab="males", xlab="females", xlim=c(70,115), ylim=c(70,115))
abline(a=0, b=1)

## Fit a model with each genotype assumed independent (ignoring pedigree structure):
fit1 <- lm(pct_D4 ~ Sex + Strain + Trt + Sex*Trt + Strain*Trt, data=model)
anova(fit1)
summary(fit1) 

## There is a sex x treatment effect, with p=0.01669,
## a strain x treatment effect, with p=2.314e-05,
## a treatment effect of p<2.2e-16,
## and a strain effect with p=4.039e-09

model1 <- aggregate(pct_D4 ~ Trt + Sex + Block + Strain, data=data, mean)

## Were males and females completed in the same block?
## Yes, for the following crosses only: (17 of 62 crosses, or 2~7% of the diallel)
## AD, AF, AH, BB, BC, CC, DA, DB, DD, DG, DH, EC, EE (some), FC, FH (some)
## GC, HF
library(beeswarm)

pdf("sex_diff_baches.pdf")
par(mfrow=c(2,2))
strain.list <- c("AD", "AF", "AH", "BB", "BC", "CC", "DA", "DB", "DD", "DG", "DH",
				 "EC", "EE", "FC", "FH", "GC", "HF")
for(i in strain.list){
	dat <- subset(data, Strain==i)
	dat <- subset(dat, (!is.na(pct_D4)))
	dat$Trt <- factor(dat$Trt)
	dat$Sex <- factor(dat$Sex)
	dat$Strain <- factor(dat$Strain)
	dat$Block <- factor(dat$Block)
	boxplot(pct_D4 ~ Strain*Trt*Sex + Block, data=dat, cex.axis=0.7, border="gray", boxwex=0.2, 
			ylim=c(70, 105), main=i)
	beeswarm(pct_D4 ~ Strain*Trt*Sex + Block, data=dat, add=TRUE, pch=16)
	dat <- subset(dat, Trt=="FLU")
	boxplot(pct_D4 ~ Strain*Sex + Block, data=dat, cex.axis=0.7, border="gray", boxwex=0.2, 
			ylim=c(70,105), main=i)
	beeswarm(pct_D4 ~ Strain*Sex + Block, data=dat, add=TRUE, pch=16)
}
dev.off()

pdf("sex_diff_baches_varylim.pdf")
par(mfrow=c(2,2))
strain.list <- c("AD", "AF", "AH", "BB", "BC", "CC", "DA", "DB", "DD", "DG", "DH",
				 "EC", "EE", "FC", "FH", "GC", "HF")
for(i in strain.list){
	dat <- subset(data, Strain==i)
	dat <- subset(dat, (!is.na(pct_D4)))
	dat$Trt <- factor(dat$Trt)
	dat$Sex <- factor(dat$Sex)
	dat$Strain <- factor(dat$Strain)
	dat$Block <- factor(dat$Block)
	boxplot(pct_D4 ~ Strain*Trt*Sex + Block, data=dat, cex.axis=0.7, border="gray", boxwex=0.2, 
			main=i)
	beeswarm(pct_D4 ~ Strain*Trt*Sex + Block, data=dat, add=TRUE, pch=16)
	dat <- subset(dat, Trt=="FLU")
	boxplot(pct_D4 ~ Strain*Sex + Block, data=dat, cex.axis=0.7, border="gray", boxwex=0.2, 
			main=i)
	beeswarm(pct_D4 ~ Strain*Sex + Block, data=dat, add=TRUE, pch=16)
}
dev.off()

pdf("sex_diff_all.pdf")
par(mfrow=c(2,1))
strain.list <- sort(unique(data$Strain))
for(i in strain.list){
	dat <- subset(data, Strain==i)
	dat <- subset(dat, (!is.na(pct_D4)))
	dat$Trt <- factor(dat$Trt)
	dat$Sex <- factor(dat$Sex)
	dat$Strain <- factor(dat$Strain)
	dat$Block <- factor(dat$Block)
	boxplot(pct_D4 ~ Block + Sex*Strain*Trt, data=dat, cex.axis=0.4, border="gray", boxwex=0.2, ylim=c(70, 115),
			main=i)
	beeswarm(pct_D4 ~ Block + Sex*Strain*Trt, data=dat, add=TRUE, pch=16)
	dat <- subset(dat, Trt=="FLU")
	boxplot(pct_D4 ~ Block + Strain*Sex, data=dat, cex.axis=0.4, border="gray", boxwex=0.2, main=i)
	beeswarm(pct_D4 ~ Block + Strain*Sex, data=dat, add=TRUE, pch=16)
}
dev.off()

pdf("sex_diff_BD_DB.pdf")
par(mfrow=c(2,1))
strain.list <- c("BB", "BD", "DB", "DD")
for(i in strain.list){
	dat <- subset(data, Strain==i)
	dat <- subset(dat, (!is.na(pct_D4)))
	dat$Trt <- factor(dat$Trt)
	dat$Sex <- factor(dat$Sex)
	dat$Strain <- factor(dat$Strain)
	dat$Block <- factor(dat$Block)
	boxplot(pct_D4 ~ Block + Sex*Strain*Trt, data=dat, cex.axis=0.4, border="gray", boxwex=0.2, ylim=c(70, 115),
			main=i)
	beeswarm(pct_D4 ~ Block + Sex*Strain*Trt, data=dat, add=TRUE, pch=16)
	dat <- subset(dat, Trt=="FLU")
	boxplot(pct_D4 ~ Block + Strain*Sex, data=dat, cex.axis=0.4, border="gray", boxwex=0.2, main=i)
	beeswarm(pct_D4 ~ Block + Strain*Sex, data=dat, add=TRUE, pch=16)
}
dev.off()

## Were reciprocals completed in the same block?
## Yes, for the following reciprocals only:

strain.list <- sort(unique(m1$Strain))
recips.test <- NULL
strain.recips <- NULL
for(x in strain.list){
	x.rev <- paste(substr(x,2,2), substr(x,1,1), sep="")
	## print(paste(x, x.rev, sep='-'))
	if(x.rev %in% strain.list & (!(x.rev==x))){
			recips.test <- c(recips.test, 1); strain.recips <- c(strain.recips, x.rev)
			}else{
			recips.test <- c(recips.test, 0); strain.recips <- c(strain.recips, NULL)}
}

recip.diff <- NULL
strain.ref <- NULL

for(i in strain.recips){
	i.rev <- paste(substr(i,2,2), substr(i,1,1), sep="")
	strain.ref <- c(strain.ref, i.rev)
	dat <- subset(m1, Strain %in% c(i, i.rev))
	recip.diff <- c(recip.diff, abs(subset(dat, Strain==i)$pct_m_minus_f) + abs(subset(dat, Strain==i.rev)$pct_m_minus_f))
}

unique.recips <- sort(unique(c(strain.ref, strain.recips)))

recip.dat <- data.frame(Strain=strain.ref, Recip=strain.recips, Sum.of.Sex.Diff=recip.diff)
recip.dattab <- recip.dat[order(-recip.diff),]
diff <- recip.dattab$Sum.of.Sex.Diff
x <- length(diff)
plot(diff ~ c(1:x), xlab="Reciprocals (ordered)", xaxt='n')
## This has twice the number of points I would want; drop every other data row

ind <- seq(1, nrow(recip.dattab), by=2)
recip.dat1 <- recip.dattab[ind, ]

diff <- recip.dat1$Sum.of.Sex.Diff
x <- length(diff)

pdf("recip_diff.pdf")
plot(diff, c(1:x), xlab="Reciprocals (ordered)", xaxt='n', ylab="abs(M-F) for AxB plus abs(M-F) for BxA", xlim=c(0,20))
text(diff, c(1:x), paste(recip.dat1$Strain, recip.dat1$Recip, sep="+"), cex=0.6, pos=4, col="red")
dev.off()

pdf("DD.pdf")
beeswarm(temp$pct_D4~temp$Trt, main="DD", ylim=c(70, 105), ylab="", xlab="white circle = female, black dot = male")
beeswarm(pct_D4~Trt, data=subset(temp, Sex=="M"), main="DD", ylim=c(70, 105), ylab="", add=TRUE, pch=16)
dev.off()
