library(plyr)
library(ggplot2)
library(gridExtra)
library(treatmentResponseDiallel)

PSq1 <- read.csv("./imputation-1000-D1/output/PSqTableStacked.csv", header = TRUE)
PSq2 <- read.csv("./imputation-1000-D2/output/PSqTableStacked.csv", header = TRUE)
PSq3 <- read.csv("./imputation-1000-D3/output/PSqTableStacked.csv", header = TRUE)
PSq4 <- read.csv("./imputation-1000-D4/output/PSqTableStacked.csv", header = TRUE)

# Note that user specified fixed effects and random effects will be 0, 
# but that the fixed effects of sex, inbreeding and sex-by-inbreeding will not.

PSq1$id <- "D1"
PSq2$id <- "D2"
PSq3$id <- "D3"
PSq4$id <- "D4"
PSq <- rbind(PSq1, PSq2, PSq3, PSq4)
PSq$id <- factor(PSq$id, levels = c("D4", 
                                    "D3",
                                    "D2", 
                                    "D1"))
PSq$X <- simplify.labels(PSq$X)
PSq <- subset(PSq, X %in% c("additive", "maternal", 
              							"inbreeding", "symmetric.epistatic", "asymmetric.epistatic",
              						  "sex-by-additive","sex-by-maternal","sex-by-inbreeding", 
              							"sex-by-symmetric.epistatic", 	"sex-by-asymmetric.epistatic"))
PSq$X <- factor(PSq$X, levels=c("additive", "maternal", 
                                "inbreeding", "symmetric.epistatic", "asymmetric.epistatic",
                                "sex-by-additive","sex-by-maternal","sex-by-inbreeding", 
                                "sex-by-symmetric.epistatic", 	"sex-by-asymmetric.epistatic"))
cbPalette <- c(	"#009E73", "#E69F00", "#56B4E9", "#0072B2",
				"#CC79A7",
				"#cb4b16",
				"#6c71c4", "#268bd2",
				"#2aa198", "#859900")

PSq$revX <- factor(PSq$X, levels = rev(levels(PSq$X)))
PSq$revid <- factor(PSq$id, levels = rev(levels(PSq$id)))
PSq$ordernum <- rev(c(1,6,2,7,3,8,4,9,5,10))
PSq$order <- factor(as.character(PSq$ordernum), levels=as.factor(c(1:10)))
PSq$color <- cbPalette[PSq$ordernum]
PSq <- PSq[with(PSq, order(revid, ordernum)),]

# ------------------------------------------------------------------------------
# Plot PSq mean and 95% credibility intervals
# ------------------------------------------------------------------------------
p1 <- ggplot(data = PSq, aes(x = order, y = mean, ymin = lower.bound, ymax = upper.bound,
                                               group = revX, colour = revX, order=revX)) +
  geom_hline(yintercept = 0, col = "grey", lwd = 1, linetype = 1) +
  geom_pointrange(stat = "identity", position = position_dodge(1), size = 0.5) + 
  coord_flip() +
  scale_color_manual(values = cbPalette[c(10:1)], name = "") +
  facet_grid(revid ~ ., switch = "y", labeller = label_wrap_gen(width = 9.5)) +
  scale_y_continuous(breaks = seq(-0.1, 1, by = 0.1), limits = c(-0.02, 0.55), expand=c(0,0)) +
  ylab("Variance Projection (CrI)") +
  theme(legend.position = "none", #panel.border = element_rect(colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        text = element_text(size = 12), 
        axis.line.x = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

# ------------------------------------------------------------------------------
# Stacked barplot - proportion of variation attributed to each inheritance class
# ------------------------------------------------------------------------------
# create a dummy variable for plotting
PSq$dummy <- 1

# ggplot doesn't like negative values for geom_bar
#if(PSq$mean < 0){PSq$mean <- 0}

p2 <- ggplot(PSq, aes(x = dummy, y = mean, fill = revX, order=revX)) +
  geom_hline(yintercept = 0, col = "grey", lwd = 1, linetype = 1) +
  geom_bar(stat = "identity") + 
  coord_flip() +
  ylab("Variance Projection (stacked)") +
  xlab("") +
  scale_fill_manual(values = cbPalette[c(10:1)], name = "") +
  scale_y_continuous(breaks = seq(-0.1, 1, by = 0.1), limits = c(-0.02, 0.55), expand=c(0,0)) +
  facet_grid(revid ~., labeller = label_wrap_gen(width = 9.5)) + 
  theme(legend.position = "none", #panel.border = element_rect(colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        text = element_text(size = 12), 
        axis.line.x = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) 

# ------------------------------------------------------------------------------
# create and extract legend
# ------------------------------------------------------------------------------
legend <- ggplot(PSq, aes(x = id, y = mean, fill = X)) +  geom_bar(stat = "identity") + 
  scale_fill_manual(values = cbPalette[c(1:10)], name = "") +
  theme(legend.position = "bottom", text = element_text(size = 13))

# function to extract legend
# source: https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

mylegend <- g_legend(legend)

# ------------------------------------------------------------------------------
# Combine plots!
# ------------------------------------------------------------------------------
pdf("VarPsAcrossTime.pdf", width=13, height=5)
grid.arrange(p1, p2, mylegend, ncol = 2, nrow = 2, 
             layout_matrix = rbind(c(1,2), c(3,3)),
             widths = c(2, 3), heights = c(2, 0.25))
dev.off()

# ------------------------------------------------------------------------------
## Varps for Mx-corrected
# ------------------------------------------------------------------------------

PSq1 <- read.csv("./imputation-1000-D1-Diplo6/output/PSqTableStacked.csv", header = TRUE)
PSq2 <- read.csv("./imputation-1000-D2-Diplo6/output/PSqTableStacked.csv", header = TRUE)
PSq3 <- read.csv("./imputation-1000-D3-Diplo6/output/PSqTableStacked.csv", header = TRUE)
PSq4 <- read.csv("./imputation-1000-D4-Diplo6/output/PSqTableStacked.csv", header = TRUE)

PSq1$id <- "D1"
PSq2$id <- "D2"
PSq3$id <- "D3"
PSq4$id <- "D4"
PSq <- rbind(PSq1, PSq2, PSq3, PSq4)
PSq$id <- factor(PSq$id, levels = c("D4", 
                                    "D3",
                                    "D2", 
                                    "D1"))
PSq$X <- simplify.labels(PSq$X)
PSq <- subset(PSq, X %in% c("additive", "maternal", 
                            "inbreeding", "symmetric.epistatic", "asymmetric.epistatic",
                            "sex-by-additive","sex-by-maternal","sex-by-inbreeding", 
                            "sex-by-symmetric.epistatic", 	"sex-by-asymmetric.epistatic"))
PSq$X <- factor(PSq$X, levels=c("additive", "maternal", 
                                "inbreeding", "symmetric.epistatic", "asymmetric.epistatic",
                                "sex-by-additive","sex-by-maternal","sex-by-inbreeding", 
                                "sex-by-symmetric.epistatic", 	"sex-by-asymmetric.epistatic"))
cbPalette <- c(	"#009E73", "#E69F00", "#56B4E9", "#0072B2",
                "#CC79A7",
                "#cb4b16",
                "#6c71c4", "#268bd2",
                "#2aa198", "#859900")

PSq$revX <- factor(PSq$X, levels = rev(levels(PSq$X)))
PSq$revid <- factor(PSq$id, levels = rev(levels(PSq$id)))
PSq$ordernum <- rev(c(1,6,2,7,3,8,4,9,5,10))
PSq$order <- factor(as.character(PSq$ordernum), levels=as.factor(c(1:10)))
PSq$color <- cbPalette[PSq$ordernum]
PSq <- PSq[with(PSq, order(revid, ordernum)),]

# ------------------------------------------------------------------------------
# Plot PSq mean and 95% credibility intervals
# ------------------------------------------------------------------------------
p1 <- ggplot(data = PSq, aes(x = order, y = mean, ymin = lower.bound, ymax = upper.bound,
                                               group = revX, colour = revX, order=revX)) +
  geom_hline(yintercept = 0, col = "grey", lwd = 1, linetype = 1) +
  geom_pointrange(stat = "identity", position = position_dodge(1), size = 0.5) + 
  coord_flip() +
  scale_color_manual(values = cbPalette[c(10:1)], name = "") +
  facet_grid(revid ~ ., switch = "y", labeller = label_wrap_gen(width = 9.5)) +
  scale_y_continuous(breaks = seq(-0.1, 1, by = 0.1), limits = c(-0.02, 0.55), expand=c(0,0)) +
  ylab("Variance Projection (CrI)") +
  theme(legend.position = "none", #panel.border = element_rect(colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        text = element_text(size = 12), 
        axis.line.x = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

# ------------------------------------------------------------------------------
# Stacked barplot - proportion of variation attributed to each inheritance class
# ------------------------------------------------------------------------------
# create a dummy variable for plotting
PSq$dummy <- 1

# ggplot doesn't like negative values for geom_bar
#if(PSq$mean < 0){PSq$mean <- 0}

p2 <- ggplot(PSq, aes(x = dummy, y = mean, fill = revX)) + #revX
  geom_hline(yintercept = 0, col = "grey", lwd = 1, linetype = 1) +
  geom_bar(stat = "identity") + 
  coord_flip() +
  ylab("Variance Projection (stacked)") +
  xlab("") +
  scale_fill_manual(values = cbPalette[c(10:1)], name = "") +
#  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(-0.1, 1)) +
  scale_y_continuous(breaks = seq(-0.1, 1, by = 0.1), limits = c(-0.02, 0.55), expand=c(0,0)) +
  facet_grid(revid ~., labeller = label_wrap_gen(width = 9.5)) + 
  theme(legend.position = "none", #panel.border = element_rect(colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        text = element_text(size = 12), 
        axis.line.x = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) 

# ------------------------------------------------------------------------------
# create and extract legend
# ------------------------------------------------------------------------------
legend <- ggplot(PSq, aes(x = id, y = mean, fill = X)) +  geom_bar(stat = "identity") + 
  scale_fill_manual(values = cbPalette[c(1:10)], name = "") +
  theme(legend.position = "bottom", text = element_text(size = 13))

# function to extract legend
# source: https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

mylegend <- g_legend(legend)

# ------------------------------------------------------------------------------
# Combine plots!
# ------------------------------------------------------------------------------
pdf("VarPsAcrossTime-Diplo6.pdf", width=13, height=5)
grid.arrange(p1, p2, mylegend, ncol = 2, nrow = 2, 
             layout_matrix = rbind(c(1,2), c(3,3)),
             widths = c(2, 3), heights = c(2, 0.25))
dev.off()
