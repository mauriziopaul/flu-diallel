library("treatmentResponseDiallel")
library("tools")

PSq <- read.csv("./D0/PSqTable.csv")
PSq$X <- simplify.labels(PSq$X)
PSq <- subset(PSq, X %in% c("additive", "maternal", 
             "inbreed.overall", "inbreeding", "symmetric.epistatic", "asymmetric.epistatic",
             "female.overall", "sex-by-additive","sex-by-maternal", "female.inbred", "sex-by-inbreeding", 
             "sex-by-symmetric.epistatic", 	"sex-by-asymmetric.epistatic", "total.explained", "Noise"))
PSq$X <- factor(PSq$X, levels=c("additive", "maternal", 
             "inbreed.overall", "inbreeding", "symmetric.epistatic", "asymmetric.epistatic",
             "female.overall", "sex-by-additive","sex-by-maternal", "female.inbred", "sex-by-inbreeding", 
             "sex-by-symmetric.epistatic", 	"sex-by-asymmetric.epistatic", "total.explained", "Noise"))

rownames(PSq) <- PSq$X

PSq <- PSq[c("additive", "maternal", 
             "inbreed.overall", "inbreeding", "symmetric.epistatic", "asymmetric.epistatic",
             "female.overall", "sex-by-additive","sex-by-maternal", "female.inbred", "sex-by-inbreeding", 
             "sex-by-symmetric.epistatic", 	"sex-by-asymmetric.epistatic", "total.explained", "Noise"),]

write.csv(PSq, "./D0/PSqTableReordered.csv")


plotMips <- function (fdir, fname, reorder = FALSE, ...)
{
    mip <- read.csv(file.path(fdir, fname))
    pdf(file.path(fdir, paste(file_path_sans_ext(fname), "_fig.pdf",
        sep = "")), width = 10, height = 6)
    par(mar = c(5.1, 13.5, 2.1, 0.1))
    if (reorder == TRUE) {
    	rownames(mip) <- mip$X
        reorder <- c("ProbFixed:Mu", #"ProbFixed:FixedEffect:1",
            "Prob:tau:RandomEffect:1", "Prob:tau:aj", "Prob:tau:motherj",
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
    mip$X[2] <- "random: block"
    rows <- nrow(mip)
    xx <- barplot(height = rev(mip$x), horiz = TRUE, names.arg = rev(mip$X),
        las = 1, xlab = "model inclusion probability (MIP) \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t",
        xlim = c(-0.01, 1.19))
    abline(v = c(0, 0.01, 0.02, 0.25, 0.5, 0.75, 0.98, 0.99,
        1), lty = 2, col = "black")
    abline(v = 0.5, lty = 1, col = "black")
    text(y = xx, x = 1.01, label = as.character(sprintf("%.3f",
        round(rev(mip$x), 3))), pos = 4)
    #try(segments(mip$x - mip$sd, rev(xx), mip$x + mip$sd, rev(xx),
    #    lwd = 1.5))
    dev.off()
}

mips <- read.csv("./D0/MIP.csv")

plotMips(fdir=savedir, fname="MIP.csv", reorder=TRUE)


plotVarps <- function (fdir, fname, reorder = FALSE, ...)
{
    psq <- read.csv(file.path(fdir, fname), row.names=1)
    psq <- psq[!grepl("FixedEffect", psq$X), ]
    psq <- psq[!grepl("RandomEffect", psq$X), ]
    if (reorder == TRUE) {
    	rownames(psq) <- psq$X
        psq <- psq[c("aj", "motherj", "BetaInbred:Av", "dominancej",
            "SymCrossjk", "ASymCrossjkDkj", "Gender:Av", "Gender:aj",
            "Gender:motherj", "BetaInbred:Gender:Av", "Gender:dominancej",
            "Gender:SymCrossjk", "Gender:ASymCrossjkDkj", "total.explained",
            "Noise"), ]
    }
    pdf(file.path(fdir, paste(file_path_sans_ext(fname), "_fig.pdf",
        sep = "")), width = 10, height = 6)
    par(mar = c(5.1, 13.5, 2.1, 12.1))
    rows <- nrow(psq)
    plot(y = c(0.5, rows + 0.5), x = c(-0.01, 1.01), yaxs = "i",
        xaxs = "i", col = "white", ylab = "", xlab = "VarP C.I.",
        yaxt = "n")
    abline(v = c(0, 0.2, 0.4, 0.6, 0.8, 1), lty = 2, col = "grey")
    abline(h = 2.5, lwd = 1.5)
    segments(x0 = psq$lower.bound, x1 = psq$upper.bound, y0 = c(rows:1),
        y1 = c(rows:1), lwd = 2)
    axis(side = 2, at = c(rows:1), labels = as.vector(simplify.labels(psq$X)),
        las = 1)
    rightlabels <- paste(as.character(sprintf("%.2f", round(100 *
        psq$mean, digits = 2))), "% (", as.character(sprintf("%.2f",
        round(100 * psq$lower.bound, digits = 2))), "%, ", as.character(sprintf("%.2f",
        round(100 * psq$upper.bound, digits = 2))), "%)", sep = "")
    axis(side = 4, at = c(rows:1), labels = rightlabels, tick = FALSE,
        las = 1)
    points(x = psq$mean, y = c(rows:1), pch = 16, col = "black",
        cex = 1.5)
    dev.off()
}

#plotVarps(fdir=savedir, fname="PSqTableReordered.csv")
