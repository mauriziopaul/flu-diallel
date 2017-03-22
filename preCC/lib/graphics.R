
A4PortraitWidth <- function(height)
{
	height/sqrt(2)
}

A4PortraitHeight <- function(width)
{
	sqrt(2)*width
}


build.mat <- function(x, y, g)
{
    ux <- sort(unique(x))
    ug <- sort(unique(g))
    mat <- matrix(nrow=length(ux), ncol=length(ug), )

    for (ix in 1:length(ux))
    {
        for (ig in 1:length(ug))
        {
            i <- which(x==ux[ix] & g==ug[ig])
            if (0==length(i)) next
            if (1<length(i))
            {
                stop("Cannot process replicate combinations: ",ux[ix],",",ug[ig],"\n")
            }
            mat[ix,ig] <- y[i]
        }
    }
    rownames(mat) <- ux
    colnames(mat) <- ug
    mat
}

plotmany <- function(formula, data=NULL, ...)
{
    form.string <- deparse(formula)
    y.form <- sub("~.*","", form.string)
    x.form <- sub("\\|.*", "", sub(".*~","", form.string))
    g.form <- sub(".*~.*\\|", "", form.string)

    y <- apply.transform(y.form, data)
    x <- apply.transform(x.form, data)
    g <- apply.transform(g.form, data)

    mat <- build.mat(x,y,g)
    thing <- matplot(as.numeric(rownames(mat)), mat, ...)
    invisible(mat)
}

my.truehist <- function (data, nbins = "Scott", h, x0 = -h/1000, breaks, prob = TRUE,
    xlim = range(breaks), ymax = max(est), col = 5, xlab = deparse(substitute(data)),
    bty = "n", add=FALSE, ...)
# truehist() from library MASS but with an "add" argument
{
    plot.truehist <- function(breaks, est, xlim, ymax, bty, xlab,
        ylab = "", density = NULL, angle = 45, col = NULL, border = NULL,
        lty = NULL, lwd = par("lwd"), ...) {
        n <- length(breaks)

        if (!add)
        {
            plot(xlim, c(0, ymax), type = "n", xlab = xlab, ylab = ylab,
                bty = bty, ...)
        }
        rect(breaks[-n], 0, breaks[-1], est, density = density,
            angle = angle, col = col, border = border, lty = lty,
            lwd = lwd)
    }
    xlab
    data <- data[is.finite(data)]
    if (missing(breaks)) {
        if (missing(h)) {
            if (is.character(nbins))
                nbins <- switch(casefold(nbins), scott = nclass.scott(data),
                  "freedman-diaconis" = , fd = nclass.FD(data))
            h <- diff(pretty(data, nbins))[1]
        }
        first <- floor((min(data) - x0)/h)
        last <- ceiling((max(data) - x0)/h)
        breaks <- x0 + h * c(first:last)
    }
    if (any(diff(breaks) <= 0))
        stop("breaks must be strictly increasing")
    if (min(data) < min(breaks) || max(data) > max(breaks))
        stop("breaks do not cover the data")
    db <- diff(breaks)
    if (!prob && sqrt(var(db)) > mean(db)/1000)
        warning("Uneven breaks with prob = FALSE will give a misleading plot")
    bin <- cut(data, breaks, include.lowest = TRUE)
    est <- tabulate(bin, length(levels(bin)))
    if (prob)
        est <- est/(diff(breaks) * length(data))
    plot.truehist(breaks, est, xlim, ymax, bty = bty, xlab = xlab,
        col = col, ...)
    invisible()
}

plot.binary.interact <- function(form, data, lty=1:10, col="black",
        ylim=c(0,1),
        pch=20, main=NULL, ylab=NULL, xlab=NULL,
        legend.xy  = NULL,
        legend.txt = NULL,
        uex=0.1,
        ci=NULL, ...)
# f( case ~ fact1 * fact2 , data)
{
    spl        <- split.formula(form)
    response   <- spl$response
    y          <- data[,response]
    predictors <- string.trim(unlist(strsplit(spl$predictors, "\\*")))
    x1    <- as.factor(data[,predictors[1]])
    x2    <- as.factor(data[,predictors[2]])

    n1 <- length(levels(x1))
    n2 <- length(levels(x2))
    if (is.null(main)) main <- paste(response, "~", predictors[1], "x", predictors[2])
    if (is.null(xlab)) xlab <- predictors[2]
    if (is.null(ylab)) ylab <- response

    pad <- 0.25
    plot(c(1-pad,n2+pad), c(0,1), type="n", axes=FALSE,
            main=main, ylab=ylab, xlab=xlab, ylim=ylim, ...)
    axis(1, labels=levels(x2), at=1:n2)
    axis(2, las=1)

    if (length(lty)<n1) lty <- rep(lty, length.out=n1)
    if (length(col)<n1) col <- rep(col, length.out=n1)
    if (length(pch)<n1) pch <- rep(pch, length.out=n1)

    se.factor <- ifelse(is.null(ci), 1, qnorm(1-(1-ci)/2))
    cat(se.factor,"\n")
    for (i1 in 1:length(levels(x1)))
    {
        a1 <- levels(x1)[i1]
        y.p  <- NULL
        y.se <- NULL
        for (a2 in levels(x2))
        {
            i <- x1==a1 & x2==a2
            i <- i & !is.na(x1) & !is.na(x2)
            p  <- mean(y[i], na.rm=TRUE)
            se <- p*(1-p)/sqrt(sum(i))
            push.back(y.p, p)
            push.back(y.se, se*se.factor)
        }
        points.errbar(1:n2, y.p, y.plus=y.se, y.minus=y.se,
                pch=pch[i1], col=col[i1], uex=uex)
        lines(1:n2, y.p, lty=lty[i1], col=col[i1])
    }
    if (is.null(legend.xy)) legend.xy <- c(1, ylim[2])
    if (is.null(legend.txt)) legend.txt <- levels(x1)

    legend(legend.xy[1], legend.xy[2],
            legend=legend.txt,
            pch=pch, lty=lty, col=col)
}


plot.table <- function(data, main="", indent=0.1, box=TRUE, line=TRUE, header=colnames(data), ...)
{
	plot(c(0,ncol(data)), c(0,nrow(data)), type="n", axes=FALSE, ylab="", xlab="", main=main)
	if (line) abline(h=nrow(data)-0.5)
	if (box) box()
	for (ic in 1:ncol(data))
	{
		text(ic-1+ncol(data)*indent, nrow(data):0, c(header[ic], data[,ic]),...)
	}
}

plot.text <- function(text, add=FALSE, mar=c(2,2,2,2), ...)
{
	lines <- unlist(strsplit(text, "\n"))
	n <- length(lines)
	oldmar <- par(mar=mar); on.exit(par(mar=oldmar))
	if (!add) plot(0:1,c(0,n), type="n", axes=FALSE, xlab="", ylab="")
	text(x=rep(0,n), y=(n:1)-0.5, adj=c(0,1), labels=lines, ...)
}

panel.errbar <- function(x,
        y=NULL,
        y.plus=0,
        y.minus=0,
        y.lower=y-y.plus,
        y.upper=y+y.minus,
        uex = 0.2,
        umbrella = uex*(range(x)%*%c(1,-1))/length(x),
        col="black",
        lwd=1,
        ...)
{
    if (!is.null(y)) panel.points(x,y,col=col,...)
    for (i in 1:length(x))
    {
        panel.lines(rep(x[i],2), c(y.lower[i], y.upper[i]), col=col, lwd=lwd)
        panel.lines(x[i] + umbrella*c(-1,1), rep(y.lower[i],2), col=col, lwd=lwd)
        panel.lines(x[i] + umbrella*c(-1,1), rep(y.upper[i],2), col=col, lwd=lwd)
    }
}


points.errbar <- function(x,
        y=NULL,
        y.plus=0,
        y.minus=0,
        y.lower=y-y.minus,
        y.upper=y+y.plus,
        uex = 0.2,
        umbrella = c(uex*(range(x)%*%c(1,-1))/length(x)),
        col="black",
        lwd=1,
        ...)
{
	col <- rep(col, length.out=length(y))
	segments(x, y.lower, x, y.upper, col=col, lwd=lwd)

	if (!is.na(umbrella))
	{
		segments(x - umbrella, y.lower, x + umbrella, y.lower, col=col, lwd=lwd)
		segments(x - umbrella, y.upper, x + umbrella, y.upper, col=col, lwd=lwd)
	}
	points(x,y,col=col,...)
}

rescaled.axis <- function(side,
        m=1,
        c=0,
        fun=function(x){c+m*x},
        inv=function(y){y/m-c},
        lim=NULL,
        unit="",
        at = NULL,
        labels=NULL,
        ...)
{
    if (is.null(at) | is.null(labels))
    {
        if (is.null(lim))
		{
	       if (1==side | 3==side)
	       {
	           lim <- par("usr")[c(1,2)]
	       }
	       else
	       {
	           lim <- par("usr")[c(3,4)]
	       }
		}
        labels <- pretty(sort(fun(lim)))
        at <- inv(labels)
    }
    axis(side, at=at, labels=paste(labels,unit,sep=""), ...)
}
