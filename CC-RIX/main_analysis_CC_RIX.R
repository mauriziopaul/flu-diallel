# 1) Install INLA (http://www.r-inla.org/download)
# 2) Move contents of Diploffect_experimental/src here
# 3) Install bagpipe backend with `R CMD install --clean --no-docs bagpipe.backend_0.34.tar.gz`
# 4) Updated versions available here: http://valdarlab.unc.edu/software/bagpipe/_build/html/bagpipe.html#r-bagpipe-backend
# 5) To run this script without error, I needed to make the following modification to the file
##   `/usr/local/lib/R/3.3/site-library/INLA/bin/mac/64bit/inla.run`, since the default /usr/lib version of liblzma.5.dylib
##   is the wrong version. Change the line:
##      export DYLD_LIBRARY_PATH="$RHOME_LIB:$DIR:/usr/lib:/opt/local/lib:$DYLD_LIBRARY_PATH"
##   to 
##      export DYLD_LIBRARY_PATH="$RHOME_LIB:$DIR:/usr/local/opt/xz/lib:/usr/lib:/opt/local/lib:$DYLD_LIBRARY_PATH"source("source.scripts.diploffect.R")
source("source.scripts.diploffect.R")
dir.create("plots", recursive=TRUE, showWarnings=FALSE)

####################################### CC-RIX data with all the weight loss time point information
ccrix.data <- read.csv("Flu-CC-RIX-data.csv", header=TRUE)

ccrix.data$SUBJECT.NAME <- ccrix.data$ID # Unique ids for each mouse
ccrix.data$eff.strain <- as.integer(ccrix.data$ccrix)

## Regressing out certain fixed and random effects; is.baric indicates location of infections: Heise Lab (0) or Baric Lab (1)
library(lme4)
lmer.frame <- model.frame(D7pct ~ D0 + is.baric + Mating + Date_Infected + ccrix + eff.strain + SUBJECT.NAME, data=ccrix.data) # keep ccrix in data frame
initial.lmer.fit <- lmer(D7pct ~ 1 + D0 + is.baric + (1|Mating) + (1|Date_Infected), data=ccrix.data)

lmer.frame$D7pct_resid <- lmer.frame$D7pct - lmer.frame$D0*initial.lmer.fit@beta[2] - lmer.frame$is.baric*initial.lmer.fit@beta[3] - ranef(initial.lmer.fit)$Date_Infected[lmer.frame$Date_Infected,]
lmer.frame$resid <- lmer.frame$D7pct - lmer.frame$D0*initial.lmer.fit@beta[2] - lmer.frame$is.baric*initial.lmer.fit@beta[3] - ranef(initial.lmer.fit)$Date_Infected[lmer.frame$Date_Infected,] - ranef(initial.lmer.fit)$Mating[lmer.frame$Mating,]

new.data <- subset(lmer.frame, select=c("D7pct_resid", "ccrix"))
new.data$D7pct_resid_scaled <- scale(new.data$D7pct_resid, center=TRUE, scale=TRUE)
mean.data.unscaled <- aggregate(D7pct_resid~ccrix, data=new.data, function(x) mean(x))
colnames(mean.data.unscaled)[2] <- "D7pct_resid"
mean.data.scaled <- aggregate(D7pct_resid_scaled~ccrix, data=new.data, function(x) mean(x))
colnames(mean.data.scaled)[2] <- "D7pct_resid_scaled"
count.data <- aggregate(D7pct_resid~ccrix, data=new.data, function(x) {NROW(x)})
colnames(count.data)[2] <- "NUM.OBS"
mean.data <- merge(merge(mean.data.unscaled, mean.data.scaled, by="ccrix"), count.data, by="ccrix")
mean.data$SUBJECT.NAME <- mean.data$ccrix
mean.data <- mean.data[,-1]
write.table(mean.data, "mean-Flu-U19-weights-20161018-03.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

### Make Kinship - add small noise on diagonal to make it non-singular
ids <- mean.data$SUBJECT.NAME[!is.na(mean.data$SUBJECT.NAME)]
K <- matrix(NA, nrow=length(ids), ncol=length(ids))
for(i in 1:(length(ids)-1)){
  K[i, i] <- 1
  these.strains <- unlist(strsplit(x=as.character(ids[i]), split="x", fixed=TRUE))
  for(j in (i+1):length(ids)){
    comparison.strains <- unlist(strsplit(x=as.character(ids[j]), split="x", fixed=TRUE))
    K[i, j] <- K[j, i] <- sum(these.strains %in% comparison.strains)/2
  }
}
K[length(ids), length(ids)] <- 1
colnames(K) <- rownames(K) <- ids
K <- K + diag(nrow(K))*0.000001

### Model using scaled means of raw data
D7.inla.diploffect.UNC27478095.model3.scaled_means <- run.diploffect.inla.through.genomecache(formula=D7pct_resid_scaled ~ 1 + locus.full, data=mean.data, K=K,
                                                                                              genomecache="./segments_happy_format/", locus="UNC27478095",
                                                                                              num.draws=100, 
                                                                                              seed=1, gamma.rate=1, impute.on="SUBJECT.NAME", weights.on="NUM.OBS",
                                                                                              do.augment.of.cache=FALSE)
D7.inla.diploffect.UNC27478095.model3.scaled_means.summary <- run.diploffect.inla.summary.stats(D7.inla.diploffect.UNC27478095.model3.scaled_means)
saveRDS(D7.inla.diploffect.UNC27478095.model3.scaled_means, "D7.inla.diploffect.UNC27478095.model3.scaled_means.RDS")
saveRDS(D7.inla.diploffect.UNC27478095.model3.scaled_means.summary, "D7.inla.diploffect.UNC27478095.model3.scaled_means.summary.RDS")
#D7.inla.diploffect.UNC27478095.model3.means_of_scaled.summary <- readRDS("D7.inla.diploffect.UNC27478095.model3.scaled_means.summary.RDS")
#D7.inla.diploffect.UNC27478095.model3.scaled_means.summary <- readRDS("D7.inla.diploffect.UNC27478095.model3.scaled_means.summary.RDS")

#pdf("plots/cc-rix_scaled_means_effect_plots.pdf", height=10, width=7.5)
pdf("plots/cc-rix_scaled_means_effect_plots.pdf", height=15, width=15)
par(mfrow=c(2,2))
plot.straineff.ci(D7.inla.diploffect.UNC27478095.model3.scaled_means.summary, flip=FALSE)
plot.diplotype.ci(D7.inla.diploffect.UNC27478095.model3.scaled_means.summary, flip=FALSE)
plot.deviation.ci(D7.inla.diploffect.UNC27478095.model3.scaled_means.summary, flip=FALSE)
plot.varexp.ci(D7.inla.diploffect.UNC27478095.model3.scaled_means.summary, add.numbers=TRUE)
dev.off()



