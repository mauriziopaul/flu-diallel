#data(FluDiData)
#data(PiximusData)
#write.csv(FluDiData, file="FluDiData.csv")
library("treatmentResponseDiallel")
data(PiximusData)
filename <- "FluDiData.csv"
args <- c("All", "D0", "Block", "NULL", "pre", "TRUE", 3000, 5, 500, "./D0/")
treatment 		<- args[1]
phenotype  		<- args[2]
random 			<- "Block"
fixed 			<- args[4]
type			<- args[5]
BS 				<- args[6]
lengthChains 	<- as.numeric(args[7])
numChains 		<- as.numeric(args[8])
burnin			<- as.numeric(args[9])
savedir 		<- args[10]
ZeroOutBeforeSample <- 0
returned <- run.tr.diallel(filename=filename, savedir=savedir, treatment=treatment, phenotype=phenotype,
	random=random, fixed=fixed, type=type, BS=BS, thin=thin, lengthChains=lengthChains,
	numChains=numChains, burnin=burnin)
plotdir <- returned[[1]]
AFD <- returned[[2]]
inner.plotter(AFD=AFD, plotdir=plotdir)
plotVarps(fdir=savedir, fname="PSqTable.csv", reorder=TRUE)