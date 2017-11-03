#' @title mcmc.add.empty: Insert empty values into MCMC for plotting purposes
#' @description Useful for plot.hpd spacing.
#' 
#' @param mcmc.obj This is the mcmc object.
#' @return Returns new mcmc object with extra unlabeled parameter
#' @examples
#' ## This example may take a couple of minutes
#' data(FluDiData)
#' data(PiximusData)
#' write.csv(FluDiData, file="FluDiData.csv")
#' filename <- "FluDiData.csv"
#' args <- c("All", "D0", "NULL", "NULL", "pre", "FALSE")
#' treatment 		<- args[1]
#' phenotype  		<- args[2]
#' random 			<- args[3]
#' fixed 			<- args[4]
#' type				<- args[5]
#' BS 				<- args[6]
#' returned <- run.tr.diallel(filename=filename, treatment=treatment, phenotype=phenotype,
#'	random=random, fixed=fixed, type=type, BS=BS)
#' plotdir <- returned[[1]]
#' TreatDi <- returned[[2]]
#' inner.plotter(TreatDi=TreatDi, plotdir=plotdir)
#' @export

mcmc.add.empty <- function(mcmc.obj, ...){
	fake.col <- rep(0, length(as.matrix(mcmc.obj[[1]][,1])))
	new.obj <- as.mcmc(cbind(as.matrix(mcmc.obj), fake.col))
	return(new.obj) # have to re-bind all chains
}
