# Example:
# ped.file <- "chr1.ped"
# cages.file <- "cages.txt"
# 
# # simple
# d <- read.table(ped.file, stringsAsFactors=FALSE)[,2:4]
# colnames(d) <- c("subject", "father", "mother")
# g <- make.pedigree.graph(d$subject, d$father, d$mother)
# plot.pedigree.graph(g, vertex.size=1)
# 
# # use cages
# cages <- as.matrix(read.table(cages.file, stringsAsFactors=FALSE)[,2:3])
# cages <- as.integer(cages[match(d$subject, cages[,1]),2])
# cages[is.na(cages)] <- 0
# 
# plot.pedigree.graph(g, xpos=cages, vertex.size=1)


make.pedigree.graph <- function(subjects, fathers, mothers,
		generation = NULL,
		directed=TRUE)
{
	require(igraph)
	subjects <- as.character(subjects)
	mothers   <- as.character(mothers)
	fathers   <- as.character(fathers)
	
	edges <- rbind(cbind(subjects, mothers), cbind(subjects, fathers))
	g <- graph.edgelist(edges, directed=directed)
	
	subjectAsEdgeIndex <- match(subjects, V(g)$name)
	edgeIndexAsSubject <- match(V(g)$name, subjects)
	V(g)$subject.index <- edgeIndexAsSubject

	nv      <- length(V(g))
	ns      <- length(subjects)
	is.root <- V(g)$name=="0"
	V(g)$is.root <- is.root
	id.root <- which(is.root)-1

	V(g)$generation <- rep(0, nv)
	if (is.null(generation))
	{
		V(g)$generation <- shortest.paths(g, id.root)
	}
	else
	{
		V(g)$generation[subjectAsEdgeIndex] <- generation
	}		
	g
}

plot.pedigree.graph <- function(g,
	xpos	   			= NULL,
	xpos.even.spacing 	= TRUE,
	vertex.color		="gray",
	edge.color  		= "darkgray",
	invisible.color		="white",
	labels				= "",
	plotit				= TRUE,
	root.elevation		= 1.2,
	tk                  = FALSE,
	axes				= TRUE,
	...)
{
	nv      <- length(V(g))
	is.root <- V(g)$name=="0"
	V(g)$is.root <- is.root
	id.root <- which(is.root)-1
	
	is.vertex.subject <- !is.na(V(g)$subject.index)
	vertex.subjects <- V(g)$subject.index[is.vertex.subject]
	
	if (!is.null(xpos))
	{
		V(g)$xpos <- rep(0.5, nv)
		V(g)$xpos[is.vertex.subject]  <- xpos[vertex.subjects]
	}
	else
	{
		V(g)$xpos <- order(V(g)$name)
	}

	V(g)$label <- ifelse(all(""==labels), "", "root")
	V(g)$label[is.vertex.subject] <- labels[vertex.subjects]

	# colors
	V(g)$color <- invisible.color
	if(1<length(vertex.color))
	{
		V(g)$color[is.vertex.subject] <- vertex.color[vertex.subjects]
	}
	else
	{
		V(g)$color[is.vertex.subject] <- vertex.color
	}
	E(g)$color  <- edge.color

	# hide annoying root vertex as much as possible
#	E(g)[ to==id.root ]$color   <- invisible.color
#	E(g)[ from==id.root ]$color <- invisible.color

	xlim <- range(V(g)$xpos)
	if (xpos.even.spacing)
	{
		xlim <- c(0, max(table(V(g)$generation)))
	}
	max.gen <- max(V(g)$generation)
	ylim <- c(0, max.gen)

	# make position matrix for graph layout, set root to middle top
	posmat <- cbind(rep(mean(xlim),nv), rep(max.gen,nv))
	posmat[is.root,] <- c(mean(xlim), max(ylim)*root.elevation)
	
	space.evenly <- function(n, w) { (1:n)*(w/(n+1)) }
	
	for (gen in sort(unique(V(g)$generation)))
	{
		igen <- gen==V(g)$generation
		if (!any(igen)) next

		wgen <- which(igen)
		if (xpos.even.spacing)
		{
			# spread subjects horizontally by their xpos
			wgen <- wgen[order(V(g)$xpos[wgen])]
			posmat[wgen,1] <- space.evenly(length(wgen), diff(xlim))
		}
		else
		{
			# position subjects exactly by their xpos
			posmat[wgen,1] <- V(g)$xpos[wgen]
		}
		
		# position vertically by generation
		posmat[wgen,2] <- (max.gen - gen)
	}
	
	if (plotit & tk)
	{
		tkplot(g, vertex.label=V(g)$label, layout=posmat, ...)
	}
	else if (plotit & !tk)
	{
		plot(g, vertex.label=V(g)$label, layout=posmat, ...)
		gen2y <- function(gen){ (max.gen - gen)/max.gen*2-1 }
		big <- 2 # big enough to cover top of vertex
#		polygonh(c(-1,1), gen2y(0.5), big, density=NA, col=invisible.color)
		rect(xleft=-1, ybottom=1, xright=gen2y(0.5), ytop=big, density=NA, col=invisible.color)
		if (axes)
		{
			drawn.gens <- setdiff(sort(unique(V(g)$generation)), 0)-1
			gen.labels <- paste("F",drawn.gens, sep="")
			axis(lty=0, at=gen2y(drawn.gens+1), labels=gen.labels, las=1, side=2)
		}
	}
	retval <- list(igraph=g, layout=posmat)
	invisible(retval)
}
