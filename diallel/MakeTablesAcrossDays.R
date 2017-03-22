library("xtable")
library("treatmentResponseDiallel")

# Make MIP table
mip1 <- read.csv("imputation-1000-D1/output/new.mip.table.csv")
mip2 <- read.csv("imputation-1000-D2/output/new.mip.table.csv")
mip3 <- read.csv("imputation-1000-D3/output/new.mip.table.csv")
mip4 <- read.csv("imputation-1000-D4/output/new.mip.table.csv")

table <- cbind(mip1[,"x"], mip2[,"x"],
				mip3[,"x"], mip4[,"x"])

names(table) <- rep("mean", 4)

mip1$X

rownames(table) <- c("Mu", "female (overall)", "inbreeding (overall)", 
	"female inbreeding (overall)", "fixed: D0", "random: Block", 
	"additive", "sex-by-additive", "maternal", "sex-by-maternal",
	"inbreeding", "sex-by-inbreeding", "symmetric epistatic",
	"sex-by-symmetric epistatic", "asymmetric epistatic", 
	"sex-by-asymmetric epistatic")

rownames.ordered <- c("Mu", "fixed: D0", "random: Block", #"random: Mx1-Diplo",
	"inbreeding (overall)", 
	"additive", "maternal", "inbreeding", "symmetric epistatic",
	"asymmetric epistatic", "female (overall)", "female inbreeding (overall)", 
	"sex-by-additive", "sex-by-maternal", "sex-by-inbreeding", 
	"sex-by-symmetric epistatic", "sex-by-asymmetric epistatic")

table <- table[rownames.ordered,]

xtab <- xtable(table, digits=4)

print(xtab)

# Make MIP table with Mx1-Diplo6
mip1 <- read.csv("imputation-1000-D1-Diplo6/output/new.mip.table.csv")
mip2 <- read.csv("imputation-1000-D2-Diplo6/output/new.mip.table.csv")
mip3 <- read.csv("imputation-1000-D3-Diplo6/output/new.mip.table.csv")
mip4 <- read.csv("imputation-1000-D4-Diplo6/output/new.mip.table.csv")

table <- cbind(mip1[,"x"], mip2[,"x"],
				mip3[,"x"], mip4[,"x"])

names(table) <- rep("mean", 4)

mip1$X
#names(table) <- rep(c("mean", "sd"), 4)
rownames(table) <- c("Mu", "female (overall)", "inbreeding (overall)", 
	"female inbreeding (overall)", "fixed: D0", "random: Block", "random: Mx1-Diplo",
	"additive", "sex-by-additive", "maternal", "sex-by-maternal",
	"inbreeding", "sex-by-inbreeding", "symmetric epistatic",
	"sex-by-symmetric epistatic", "asymmetric epistatic", 
	"sex-by-asymmetric epistatic")

rownames.ordered <- c("Mu", "fixed: D0", "random: Block", "random: Mx1-Diplo",
	"inbreeding (overall)", 
	"additive", "maternal", "inbreeding", "symmetric epistatic",
	"asymmetric epistatic", "female (overall)", "female inbreeding (overall)", 
	"sex-by-additive", "sex-by-maternal", "sex-by-inbreeding", 
	"sex-by-symmetric epistatic", "sex-by-asymmetric epistatic")

table <- table[rownames.ordered,]

xtab <- xtable(table, digits=4)

print(xtab)

# Make TreVarP table
varp1 <- read.csv("imputation-1000-D1/output/PSqTableStacked.csv")
varp2 <- read.csv("imputation-1000-D2/output/PSqTableStacked.csv")
varp3 <- read.csv("imputation-1000-D3/output/PSqTableStacked.csv")
varp4 <- read.csv("imputation-1000-D4/output/PSqTableStacked.csv")


table <- cbind(varp1[,c("mean", "lower.bound", "upper.bound")],
			varp2[,c("mean", "lower.bound", "upper.bound")],
			varp3[,c("mean", "lower.bound", "upper.bound")],
			varp4[,c("mean", "lower.bound", "upper.bound")])

names(table) <- rep(c("mean", "0.025", "0.975"), 4)
rownames(table) <- varp1$X

table.sub <- table[!(rownames(table) %in% c("FixedEffect:1", "RandomEffect:1")), ]
rownames(table.sub) <- removeDots(simplify.labels(rownames(table.sub)))

rownames.ordered <- c('female overall','inbreed overall','female inbred',
	'additive','maternal','inbreeding','symmetric epistatic','asymmetric epistatic',
	'sex-by-additive','sex-by-maternal','sex-by-inbreeding','sex-by-symmetric epistatic','sex-by-asymmetric epistatic',
	'total explained','Noise')

table.sub <- table.sub[rownames.ordered,]

rownames(table.sub) <- c('female (overall)','inbreed (overall)','female inbreeding',
	'additive','maternal','inbreeding','symmetric epistatic','asymmetric epistatic',
	'sex-by-additive','sex-by-maternal','sex-by-inbreeding','sex-by-symm. epis.','sex-by-asymm. epis.',
	'total explained','noise')

xtab <- xtable(table.sub, digits=4)

print(xtab)


# Make TreVarP table, after accounting for Mx1
varp1 <- read.csv("imputation-1000-D1-Diplo6/output/PSqTableStacked.csv")
varp2 <- read.csv("imputation-1000-D2-Diplo6/output/PSqTableStacked.csv")
varp3 <- read.csv("imputation-1000-D3-Diplo6/output/PSqTableStacked.csv")
varp4 <- read.csv("imputation-1000-D4-Diplo6/output/PSqTableStacked.csv")


table <- cbind(varp1[,c("mean", "lower.bound", "upper.bound")],
			varp2[,c("mean", "lower.bound", "upper.bound")],
			varp3[,c("mean", "lower.bound", "upper.bound")],
			varp4[,c("mean", "lower.bound", "upper.bound")])

names(table) <- rep(c("mean", "0.025", "0.975"), 4)
rownames(table) <- varp1$X

table.sub <- table[!(rownames(table) %in% c("FixedEffect:1", "RandomEffect:1")), ]
rownames(table.sub) <- removeDots(simplify.labels(rownames(table.sub)))

rownames.ordered <- c('female overall','inbreed overall','female inbred',
	'additive','maternal','inbreeding','symmetric epistatic','asymmetric epistatic',
	'sex-by-additive','sex-by-maternal','sex-by-inbreeding','sex-by-symmetric epistatic','sex-by-asymmetric epistatic',
	'total explained','Noise')

table.sub <- table.sub[rownames.ordered,]

rownames(table.sub) <- c('female (overall)','inbreed (overall)','female inbreeding',
	'additive','maternal','inbreeding','symmetric epistatic','asymmetric epistatic',
	'sex-by-additive','sex-by-maternal','sex-by-inbreeding','sex-by-symm. epis.','sex-by-asymm. epis.',
	'total explained','noise')

xtab <- xtable(table.sub, digits=4)

print(xtab)

