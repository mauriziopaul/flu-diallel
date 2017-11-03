##---------------------------------------------------------------------------------------------------------------------
## Title: CC Tools Functions for Data Wrangling with CC and CC-RIX Related Data 
## Author: Paul L. Maurizio
## Institution: University of North Carolina at Chapel Hill
## Program: Bioinformatics and Computational Biology Ph.D. Curriculum
## Advisors: Will Valdar, Mark Heise
## Date Created: 2016-08-08
## Date Updated: 2016-12-19
##---------------------------------------------------------------------------------------------------------------------

#' plmCCTools: A package for wrangling CC and CC-RIX related data.
#' The plmCCTools package provides the important function:
#' ccNameConverter
#' 
#' @docType package

#' @name ccNameConverter
NULL

#' @import cmdline
NULL

#' @import configfile
NULL

#' @import methods
NULL

#' @import tools
NULL

#' @import stringr
NULL

#' @import devtools
NULL

#' @import gridGraphics
NULL

#' @section require namespaces:

requireNamespace("cmdline", quietly=TRUE)
requireNamespace("configfile", quietly=TRUE)
requireNamespace("methods", quietly=TRUE)
requireNamespace("tools", quietly=TRUE)
requireNamespace("stringr", quietly=TRUE)
requireNamespace("devtools", quietly=TRUE)

#' @section treatmentResponseDiallel functions:

#' @title colSwitch: Define standard colors for CC founders
#' @description insert something here.
#' 
#' @param arg the argument, in the form of "AA", "BB", etc.
#' @param ... additional arguments
#' @return Defines standard CC founder colors in hex-code
#' @examples
#' col.switch("A", "B", "C", "DD", "EE", "FF")
#' @export
colSwitch <- function(arg, ...){
 	founder.colors <- c("#F0F000", "#808080", "#F08080", "#1010F0", "#00A0F0", "#00A000", "#F00000", "#9000E0")
 	ret <- NULL
 	for(i in 1:length(arg)){
 		ret[i] <- switch(as.character(arg[i]),
 				AA=founder.colors[1],
 				BB=founder.colors[2],
 				CC=founder.colors[3],
 				DD=founder.colors[4],
 				EE=founder.colors[5],
 				FF=founder.colors[6],
 				GG=founder.colors[7],
 				HH=founder.colors[8],
 				A=founder.colors[1],
 				B=founder.colors[2],
 				C=founder.colors[3],
 				D=founder.colors[4],
 				E=founder.colors[5],
 				F=founder.colors[6],
 				G=founder.colors[7],
 				H=founder.colors[8],
 				"black")
 	}
 	return(ret)
}

lettersToStrains <- function(arg, ...){
 	founder.names <- c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB")
 	ret <- NULL
 	for(i in 1:length(arg)){
 		ret[i] <- switch(as.character(arg[i]),
 				AA=founder.names[1],
 				BB=founder.names[2],
 				CC=founder.names[3],
 				DD=founder.names[4],
 				EE=founder.names[5],
 				FF=founder.names[6],
 				GG=founder.names[7],
 				HH=founder.names[8],
 				A=founder.names[1],
 				B=founder.names[2],
 				C=founder.names[3],
 				D=founder.names[4],
 				E=founder.names[5],
 				F=founder.names[6],
 				G=founder.names[7],
 				H=founder.names[8],
 				"unknown")
 	}
 	return(ret)
}

#' @title lettersToNumbers
#' @description This converts strain name letters to numbers
#' 
#' @param arg the argument, in the form of "A", "B", ... "H", or vector of these.
#' @param ... additional arguments
#' @return Defines standard CC founder numbers from single letters.
#' @examples
#' lettersToNumbers(c("A", "B", "C", "D", "E", "F"))
#' @export
lettersToNumbers <- function(arg, ...){
  ret <- NULL
  for(i in 1:length(arg)){
    ret[i] <- switch(as.character(arg[i]),
                     A=1, B=2, C=3, D=4,
                     E=5, F=6, G=7, H=8)
  }
  return(ret)
}


#' @title convertCrossNames
#' @description This converts CC cross names name letters to numbers
#' 
#' @param arg the argument, in the form of "8042", "8002", etc., or vector of these.
#' @param ... additional arguments
#' @return Returns the CC--- version of the numeric strain names in the vector. 
#' @examples
#' convertCrossNames(c("8042", "8002"))
#' @export
convertCrossNames <- function(crossnames, ...){
	ret <- rep(NA, length(crossnames))
	for(i in 1:length(crossnames)){
		ret[i] <- switch(as.character(crossnames[i]),
					"18042"="CC042",
					"8042"="CC042",
					"8002"="CC032",
					"8004"="CC024",
					"8005"="CC012",
					"8008"="CC025",
					"8010"="CC013",
					"8016"="CC028",
					"8018"="CC010",
					"18018"="CC010",
					"8021"="CC043",
					"8024"="CC016",
					"8026"="CC026",
					"8027"="CC027",
					"8031"="CC031",
					"8033"="CC033",
					"8034"="CC030",
					"8036"="CC008",
					"8043"="CC023",
					"18045"="CC045",
					"8045"="CC045",
					"8046"="CC022",
					"8048"="CC061",
					"8049"="CC038",
					"8050"="CC054",
					"8052"="CC052",
					"8054"="CC020",
					"8056"="CC056",
					"16012"="CC059",
					"16034"="CC070",
					"16072"="CC037",
					"16188"="CC004",
					"16211"="CC005",
					"16296"="CC049",
					"16441"="CC041",
					"16513"="CC019",
					"16521"="CC072",
					"16557"="CC040",
					"16680"="CC055",
					"16750"="CC006",
					"16768"="CC068",
					"16785"="CC071",
					"16912"="CC051",
					"13067"="CC003",
					"13140"="CC001",
					"13421"="CC007",
					"1515"="CC073",
					"15155"="CC039",
					"15156"="CC002",
					"1566"="CC021",
					"3015"="CC074",
					"3032"="CC017",
					"3154"="CC015",
					"3252"="CC011",
					"3260"="CC075",
					"3393"="CC057",
					"3415"="CC014",
					"3460"="CC060",
					"3564"="CC029",
					"3609"="CC018",
					"4410"="CC044",
					"477"="CC036",
					"5035"="CC035",
					"5080"="CC034",
					"15119"="CC065",
					"5119"="CC065",
					"5248"="CC048",
					"5306"="CC062",
					"5343"="CC076",
					"5346"="CC046",
					"5358"="CC058",
					"5391"="CC063",
					"5489"="CC009",
					"15489"="CC009",
					"559"="OR559",
					"5612"="CC047",
					"773"="CC053",
					"867"="CC050",
					NA)
	}
	return(ret)
}

#' @title backConvertCrossNames
#' @description This converts CC cross names
#' 
#' @param arg the argument, in the form of "CC001", "CC002", etc., or vector of these.
#' @param ... additional arguments
#' @return Returns the number version of the CC--- strain names in the vector. 
#' @examples
#' backConvertCrossNames(c("CC002", "CC001"))
#' @export
backConvertCrossNames <- function(crossnames, ...){
  ret <- rep(NA, length(crossnames))
  for(i in 1:length(crossnames)){
    ret[i] <- switch(as.character(crossnames[i]),
#                     "CC042"="18042",
                     "CC042"="8042",
                     "CC032"="8002",
                     "CC024"="8004",
                     "CC012"="8005",
                     "CC025"="8008",
                     "CC013"="8010",
                     "CC028"="8016",
                     "CC010"="8018",
#                     "CC010"="18018",
                     "CC043"="8021",
                     "CC016"="8024",
                     "CC026"="8026",
                     "CC027"="8027",
                     "CC031"="8031",
                     "CC033"="8033",
                     "CC030"="8034",
                     "CC008"="8036",
                     "CC023"="8043",
                     "CC045"="8045",
                     "CC022"="8046",
                     "CC061"="8048",
                     "CC038"="8049",
                     "CC054"="8050",
                     "CC052"="8052",
                     "CC020"="8054",
                     "CC056"="8056",
                     "CC059"="16012",
                     "CC070"="16034",
                     "CC037"="16072",
                     "CC004"="16188",
                     "CC005"="16211",
                     "CC049"="16296",
                     "CC041"="16441",
                     "CC019"="16513",
                     "CC072"="16521",
                     "CC040"="16557",
                     "CC055"="16680",
                     "CC006"="16750",
                     "CC068"="16768",
                     "CC071"="16785",
                     "CC051"="16912",
                     "CC003"="13067",
                     "CC001"="13140",
                     "CC007"="13421",
                     "CC073"="1515",
                     "CC039"="15155",
                     "CC002"="15156",
                     "CC021"="1566",
                     "CC074"="3015",
                     "CC017"="3032",
                     "CC015"="3154",
                     "CC011"="3252",
                     "CC075"="3260",
                     "CC057"="3393",
                     "CC014"="3415",
                     "CC060"="3460",
                     "CC029"="3564",
                     "CC018"="3609",
                     "CC044"="4410",
                     "CC036"="477",
                     "CC035"="5035",
                     "CC034"="5080",
                     "CC065"="15119",
                     "CC065"="5119",
                     "CC048"="5248",
                     "CC062"="5306",
                     "CC076"="5343",
                     "CC046"="5346",
                     "CC058"="5358",
                     "CC063"="5391",
                     "CC009"="5489",
                     "CC009"="15489",
                     "OR559"="559",
                     "CC047"="5612",
                     "CC053"="773",
                     "CC050"="867",
                     NA)
  }
  return(ret)
}

#' @title underNumAlpha
#' @description This adds underscores between alphabet and numeric characters
#' 
#' @param arg the argument which has numeric and character strings next to each ther
#' or a vector of character strings
#' @param ... additional arguments
#' @return Returns elements with underscores between letters and numbers
#' @examples
#' underNumAlpha(c("135RL","136R"))
#' @export
underNumAlpha <- function(ids, ...){
	ret <- rep(NA, length(ids))
		for(i in 1:length(ids)){
			num <- str_extract(ids[i], "^[0-9]+")
			alpha <- str_extract(ids[i], "[A-Z]+$")
			ret[i] <- paste(num, alpha, sep="_")
		}
	return(ret)
}

#' @title removeDots
#' @description This removes dots (and underscores) from a string, replacing with a space.
#' 
#' @param x the dot-containing string
#' or a vector of character strings
#' @param ... additional arguments
#' @return Returns string without dot/underscore.
#' @examples
#' removeDots("This_is_a.test.string.")
#' @export
removeDots <- function(x, ...){
  if(!(is.character(x))){
  	cat("ERROR:", as.character(x), "is not a string! \n")
		return(cat(NULL))
	}
	new.x <- gsub(".", " ", x, fixed=TRUE)
	new.x <- gsub("_", " ", new.x, fixed=TRUE)
	return(new.x)
}

#' @title docBuildInstall
#' @description This documents, builds, and installs an R package based on its directory
#' 
#' @param dir the directory of the R package to be documented, built, installed
#' @param ... additional arguments
#' @return Returns string without dot/underscore.
#' @examples
#' # Not Run
#' @export
docBuildInstall  <- function(dir, ...){
	document(dir, roclets=c('rd', 'namespace'))
	build(dir)
	install(dir, quick=TRUE, upgrade_dependencies=FALSE)
}


parseDamSire <- function(data.frame, cross.colname){
	data.frame$dam <- sapply(data.frame[, as.character(cross.colname)], 
		FUN=function(x){gsub("x.*", "", x)})
	data.frame$sire <- sapply(data.frame[, as.character(cross.colname)], 
		FUN=function(x){gsub(".*x", "", x)})
	return(data.frame)
}

#' @title plotCheckerboard
#' @description This plots checkerboard-style heatmaps of mean phenotype values in the diallel 
#' @param summary.dat summary data frame, summarized by strain-pair
#' @param phen phenotype column name
#' @param aggregate.col data frame column to aggregate by 
#' @param fun usually supplied as "mean", but could also be "sd" or some similar function of interest
#' @param founder.ids usually letters A through ...
#' @param founder.names descriptive names of founders, for axis labels
#' @param founder.colors colors for founder lines, for axis plotting
#' @param match.by choose founder.letters or founder.names, based on what is in input data frame
#' @param scale.palette specify the palette to be used for color scale in the heatmap
#' @param scale.range define the range of the data
#' @param plot.palette specify the palette for the scale to be plotted, which may differ from the range of the palette of the plot
#' @param scale.tick specify whether to highlight a mark on the axis with an additional tick 
#' @param which.tick specify which axis tick to highlight
#' @param noscale boolean specifying whether to hide the scale
#' @param scale.las specify the side of the scale axis using las
#' @param axis.srt specify the angle to tilt the scale axis labels
#' @param axis.y.adj specify a vertical adjustment for the scale axis labels
#' @param scale.cex specify how big the scale text is using standard cex adjustments 
#' @param ... additional arguments
#' @return Returns checkerboard plot
#' @examples
#' # Not Run
#' @export

plotCheckerboard <- function(summary.dat=NULL, phen=NULL, aggregate.col="Strain", fun="mean",
                             founder.ids=LETTERS[1:8], main="", xlab="",
                             founder.names=c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB"),
                             founder.letters=LETTERS[1:8], 
                             founder.colors=c("#F0F000", "#808080", "#F08080", "#1010F0", 
                                                 "#00A0F0", "#00A000", "#F00000", "#9000E0"),
                             match.by=c("founder.letters","founder.names")[1],
                             scale.palette=colorRampPalette(c("white", "black"))(n = 1000),
                             scale.range=range(summary.dat, na.rm=TRUE),
                             plot.palette=scale.palette, scale.tick=FALSE, which.tick=5,
                             noscale=c(FALSE, TRUE)[1], scale.las=1, axis.srt=45, axis.y.adj=0.1,
                             scale.cex=1, scale.gradient.cex=1,
                             ...){
  fdim <- length(founder.ids)
  if(match.by=="founder.letters"){
    matchvec <- expand.grid(founder.letters, founder.letters)
    matchvec$paste <- paste0(matchvec$Var1, matchvec$Var2)
  }
  if(match.by=="founder.names"){
    matchvec <- expand.grid(founder.names, founder.names)
    matchvec$paste <- paste(matchvec$Var1, matchvec$Var2, sep="x")
  }
  matchmat <- matrix(as.character(matchvec$paste), nrow=fdim, ncol=fdim)
  data.mat <- matrix(NA, nrow=fdim, ncol=fdim)
  colnames(data.mat) <- rownames(data.mat) <- founder.names
  
  for(j in 1:fdim){
    for(k in 1:fdim){
      thiscross <- as.character(matchmat[j,k])
      try(data.mat[j,k] <- as.numeric(summary.dat[summary.dat[,aggregate.col]==thiscross,phen]),
          silent=TRUE)
    }
  }
  
  heatfun <- function(){
    heatmap(t(data.mat), scale="none", Colv=NA, Rowv=NA,
            revC=TRUE, symm=TRUE, ColSideColors = founder.colors, margins=c(6,6),
            RowSideColors = founder.colors, xlab=expression(bold("sire")), 
            ylab=expression(bold("dam")), col=plot.palette, main=main)
  }
  grid.newpage()
  pushViewport(viewport(y=1.1, x=0.5, width=0.8, just="top"))
  grid.echo(heatfun, newpage=FALSE)
  upViewport()
  
  scalefun <- function(){
    plot(c(1:length(scale.palette)), rep(0.92,length(scale.palette)), 
         pch="|", col=scale.palette, cex=scale.gradient.cex, ylab="", xlab=xlab,
         yaxt="n", frame.plot=FALSE, ylim=c(0.85,1.4), xaxt="n",
         adj=0, yaxs='i')
    labs <- signif(as.numeric(pretty(scale.range)), digits=3)
    if(TRUE==scale.tick){
      axis(side=3, at=seq(from=1, to=length(scale.palette), 
                          length.out=length(labs))[which.tick],
           labels="", pos=0.9)
    }
    axis(side=1, at=seq(from=1, to=length(scale.palette), length.out=length(labs)),
         labels=FALSE, pos=0.9, las=scale.las)
    text(x=seq(from=1, to=length(scale.palette), length.out=length(labs)), labels=labs,
         y=par()$usr[3]-axis.y.adj*(par()$usr[4]-par()$usr[3]),
         srt=axis.srt, adj=1, xpd=TRUE, cex=scale.cex)
  }
  if(noscale==FALSE){
    pushViewport(viewport(y=0, height=0.4, width=0.85, x=0.44, just="bottom"))
    grid.echo(scalefun, newpage=FALSE)
    upViewport()
  }
}

#' @title rescale.palette
#' @description This rescales the palette based on where the plot range falls within the data range
#' @param plot.range the range of values for plotting (i.e., range of strain means)
#' @param data.range the range of values for the color ramp scale (i.e., range of strain means across all plots)
#' @param start.pal the starting palette
#' @param ... additional arguments
#' @return Returns rescaled palette, with same length of colors as original scale
#' @examples
#' # Not Run
#' @export

rescale.palette <- function(plot.range, data.range, color.len=1000,
                            start.pal=colorRampPalette(c("white", "black"))(n = color.len),
                            ...){
  center <- -data.range[1]
  scale <- (length(start.pal)-1)/(data.range[2]-data.range[1])
  len <- length(start.pal)-1
  rescaled.startcol <- round((plot.range[1]+center)*scale + 1)
  rescaled.endcol <- round((plot.range[2]+center)*scale + 1)
  rescaled <- colorRampPalette(c(start.pal[rescaled.startcol], start.pal[rescaled.endcol]))(n = color.len)
  return(rescaled)
}