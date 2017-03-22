#!/usr/local/bin/Rscript
library(cmdline)
library(configfile)
library(treatmentResponseDiallel)
library(tools)
library(data.table)

#confile <- "settings.config.D1"
#confile <- "settings.config.D2"
#confile <- "settings.config.D3"
#confile <- "settings.config.D4"
#confile <- "settings.config.D1-Diplo6"
#confile <- "settings.config.D2-Diplo6"
#confile <- "settings.config.D3-Diplo6"
#confile <- "settings.config.D4-Diplo6"
confile <- cmdline.string("configfile")
config <- read.configfile(confile)
source("loadConfig.R")

if(batch=="NULL"){batch <- NULL}
if(fixed=="NULL"){fixed <- NULL}
if(random=="NULL"){random <- NULL}

random <- c(random,batch)
if(length(random)>1){random <- as.list(random)}
if(length(fixed)>1){fixed <- as.list(fixed)}

# stacked.posterior <- read.delim2(file.path(savedir, "output/stackedPostTreatDiMIMPs.csv"), 
# 	check.names=FALSE, sep=";", dec=",")

# stacked.posterior <- as.mcmc(stacked.posterior)
# VarPs <- summary(stacked.posterior)
# VarPs1 <- VarPs[[1]]
# VarPs2 <- VarPs[[2]]
# colnames(VarPs2)[c(3,1,5)] <- c("mean", "lower.bound", "upper.bound") ## mean should be median
# #write.csv(VarPs2, file=file.path(savedir, "output/PSqTableStacked.csv"), row.names=TRUE)
# #reorder
# if(length(random)<2){
# 	VarPs2 <- VarPs2[c(#'FixedEffect:1','RandomEffect:1',
# 					'aj','motherj','BetaInbred:Av','dominancej','SymCrossjk','ASymCrossjkDkj',
# 					'Gender:Av',
# 					'Gender:aj','Gender:motherj',
# 					'BetaInbred:Gender:Av','Gender:dominancej','Gender:SymCrossjk',
# 					'Gender:ASymCrossjkDkj','total.explained','Noise'),]
# }else{
# 	VarPs2 <- VarPs2[c(#'FixedEffect:1','RandomEffect:1','RandomEffect:2',
# 					'aj','motherj','BetaInbred:Av','dominancej','SymCrossjk','ASymCrossjkDkj',
# 					'Gender:Av',
# 					'Gender:aj','Gender:motherj',
# 					'BetaInbred:Gender:Av','Gender:dominancej','Gender:SymCrossjk',
# 					'Gender:ASymCrossjkDkj','total.explained','Noise'),]
# } # This does not account for > 2 random effects

# fwrite(as.data.frame(VarPs2), file=file.path(savedir, "output/PSqTableStacked.csv"), row.names=TRUE)
# try(plotVarps(fdir=savedir, fname="output/PSqTableStacked.csv"))

#system("open ./imputation-1000-D*/output/PSqTableStacked_fig.pdf")

##------------
## PLOT MIPs
##------------

plotMips <- function (fdir, fname, reorder = FALSE, ...)
{
    mip <- read.csv(file.path(fdir, fname))
    pdf(file.path(fdir, paste(file_path_sans_ext(fname), "_fig.pdf",
        sep = "")), width = 10, height = 6)
    par(mar = c(5.1, 13.5, 2.1, 0.1))
    if (reorder == TRUE) {
    	rownames(mip) <- mip$X
        reorder <- c("ProbFixed:Mu", "ProbFixed:FixedEffect:1",
            "Prob:tau:RandomEffect:1", "Prob:tau:RandomEffect:2", 
            "Prob:tau:aj", "Prob:tau:motherj",
            "ProbFixed:BetaInbred:Av", "Prob:tau:dominancej",
            "Prob:tau:SymCrossjk", "Prob:tau:ASymCrossjkDkj",
            "ProbFixed:Gender:Av", "Prob:tau:Gender:aj", "Prob:tau:Gender:motherj",
            "ProbFixed:BetaInbred:Gender:Av", "Prob:tau:Gender:dominancej",
            "Prob:tau:Gender:SymCrossjk", "Prob:tau:Gender:ASymCrossjkDkj")
        rownames(mip) <- mip$X
        mip <- mip[reorder, ]
    }
    mip$X <- simplify.labels(mip$X)
    mip$X[1] <- "mu"
    mip$X[2] <- "fixed: D0"
    mip$X[3] <- "random: block"
    mip$X[4] <- "random: Mx1 diplo"
    rows <- nrow(mip)
    xx <- barplot(height = rev(mip$x), horiz = TRUE, names.arg = rev(mip$X),
        las = 1, xlab = "model inclusion probability (MIP) \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t",
        xlim = c(-0.01, 1.19))
    abline(v = c(0, 0.01, 0.02, 0.25, 0.5, 0.75, 0.98, 0.99,
        1), lty = 2, col = "black")
    abline(v = 0.5, lty = 1, col = "black")
    text(y = xx, x = 1.01, label = as.character(sprintf("%.3f",
        round(rev(mip$x), 3))), pos = 4)
    try(segments(mip$x - mip$sd, rev(xx), mip$x + mip$sd, rev(xx),
        lwd = 1.5))
#    try(points(mip$x, rev(xx), pch=1, lwd=1.5))
    dev.off()
}
mip.summary <- "new.mip.table.csv"
#rownames(mip.summary)
try(plotMips(fdir=file.path(savedir, "output"), fname=mip.summary, reorder=TRUE))

cat("Finished plotting.\n")
