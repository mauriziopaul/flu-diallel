#library(rjags)



reformat.data.for.aduse <- function(data,
        components=c("a", "d", "u", "s", "e", "sex"))
{
    required.cols <- c(
            "phenotype",
            "is.female",
            "mother.strain.name",
            "father.strain.name")
    if (!all(required.cols %in% colnames(data)))
    {
        stop("Missing columns ", paste(collapse=", ", setdiff(required.cols, colnames(data))))
    }

    strain.name   <- sort(unique(c(as.character(data$mother.strain.name), as.character(data$father.strain.name))))
    mother.strain <- as.integer(factor(data$mother.strain.name, levels=strain.name))
    father.strain <- as.integer(factor(data$father.strain.name, levels=strain.name))

    parentpair.name  <- paste(data$mother.strain.name, sep="_x_", data$father.strain.name)
    parentpair.id    <- as.integer(as.factor(parentpair.name))

    strainpair.name <- ifelse(mother.strain < father.strain,
                paste(data$mother.strain.name, sep="_x_", data$father.strain.name),
                paste(data$father.strain.name, sep="_x_", data$mother.strain.name))
    strainpair.id <- as.integer(as.factor(strainpair.name))

    parentpair.to.strainpair <- integer(unique(parentpair.id))
    parentpair.to.strainpair[parentpair.id] <- strainpair.id

    augmented.data <- cbind(data, data.frame(
            mother.strain=mother.strain,
            father.strain=father.strain,
            parentpair.name=parentpair.name,
            parentpair.id=parentpair.id,
            strainpair.name=strainpair.name,
            strainpair.id=strainpair.id))

    jags.data <- list(
            n = nrow(data),
            phenotype = data$phenotype)

    if (any(components %in% c("sex")))
    {
        jags.data$is.female <- as.integer(data$is.female)
    }
    if (any(components %in% c("a", "d", "u", "s")))
    {
        jags.data$num.strains   <- length(strain.name)
        jags.data$mother.strain <- mother.strain
        jags.data$father.strain <- father.strain
    }
    if (any(components %in% c("e")))
    {
        jags.data$num.strainpairs           <- length(unique(strainpair.id))
        jags.data$parentpair.to.strainpair  <- parentpair.to.strainpair

        jags.data$num.parentpairs <- length(unique(parentpair.id))
        jags.data$parentpair.id             <- parentpair.id
    }
    list( data = augmented.data, jags.data = jags.data )
}

gelman.converged <- function(x, min.psrfq=1.1)
{
    max(gelman.diag(x)$psrf[,2]) < min.psrfq
}


plot.ci <- function(midvals, narrow.intervals, wide.intervals,
    names=1:length(midvals),
    pos = 0,
    add=FALSE,
    xlab="Estimate",
    xlab.line=0,
    ylab="",
    yaxis=TRUE,
    name.margin=6,
    name.line=4,
    pch.midvals=19,
    col="black",
    col.midvals=col,
    type="p",
    log="",
    xrev=F,
    xlim = NULL,
    ...)
{
  nvals <- length(midvals)
  col.midvals <- rep(col.midvals, length.out=nvals)
  y.pos <- (1:nvals)-0.5
  if (!add)
  {
    lim <- range(c(wide.intervals,narrow.intervals,midvals), na.rm=TRUE)
    if (log=="x"){
      lim[lim<=0]=0.0000001
    }
    else {
      lim <- c(-1,1) * diff(lim)*0.1 + lim
    }
    if (is.null(xlim)) {
      xlim = range(lim)
    }
    if (xrev) xlim = rev(xlim)
    plot(lim, c(0,nvals+0.5), type="n", axes=FALSE,
         ylab=ylab, xlab="", log=log, xlim=xlim, ...)
    title(xlab=xlab)
    axis(1)
    axis(3, line=-1)
    if (yaxis)
    {
      ## for different plot, you probably need to change
      ## the pos parameter here
      axis(2, at=y.pos, labels=rev(names), las=1,
           lty=0, hadj=0, line=name.line, outer=F,
           pos=pos)
    }
  }
  if ("p"==type)
  {
    for (i in 1:nvals)
    {
      pos <- nvals-i + 0.5
      lines(wide.intervals[i,], rep(pos,2))
      lines(narrow.intervals[i,], rep(pos,2), lwd=3)
      points(midvals[i], pos, pch=pch.midvals, col=col.midvals[i])
    }
  }
  invisible(rev(y.pos))
}

plot.corr <- function(mat,
        range=c(-1,1),
        border.lwd=0,
        num.colors=200,
        main="correlations",
        line.col="white",
        name.margin=6,
        name.line=4,
        names=colnames(mat),
        darkness.extent=1,
        positive.increasing=TRUE,
        show.range=TRUE,
        show.sign=TRUE)
{
    num.colors <- 2*ceiling(num.colors/2)
    if (TRUE)
    {
        darkness <- c( ((num.colors/2):1)/(num.colors/2), 0, (1:(num.colors/2))/(num.colors/2))
        darkness <- darkness*darkness.extent
        colors <- gray(1-darkness)
    }
    else
    {
        colors <- rainbow(num.colors+1, start=4/6, end=0)
    }
    limits <- 0:num.colors * (diff(range)/num.colors) + min(range)
    n <- nrow(mat)
    mar <- c(5, name.margin, 4, 2)+0.1
    oldmar <- par(mar=mar); on.exit(par(mar=oldmar))
    y.offset <- 0.5

    plot(c(0,n), c(0,n+y.offset), axes=FALSE, type="n", ylab="", xlab="", main=main)
    image(0:n, c(0:n)+y.offset, mat[,n:1], zlim=range, col=colors, add=TRUE)

    if (show.sign)
    {
        x0=matrix(nrow=n, ncol=n)
        y0=matrix(nrow=n, ncol=n)
        x1=matrix(nrow=n, ncol=n)
        y1=matrix(nrow=n, ncol=n)
        for (j in 1:n)
        {
            for (k in 1:n)
            {
                x <- mat[j,k]
                x.start <- j-1
                x.end   <- j
                y.start <- y.offset + n-k
                y.end   <- y.offset + n-k+1

                x0[j,k]=x.start
                x1[j,k]=x.end
                if ((x>0 & positive.increasing) | (x<0 & !positive.increasing))
                {
                    y0[j,k]=y.start
                    y1[j,k]=y.end
                }
                else if ((x<0 & positive.increasing)|(x>0 & !positive.increasing))
                {
                    y0[j,k]=y.end
                    y1[j,k]=y.start
                }
            }
        }
        segments(c(x0), c(y0), c(x1), c(y1), col=line.col)
    }
    if (0<border.lwd)
    {
        abline(h=(0:n)+y.offset, lwd=border.lwd, col="white")
        abline(v=0:n, lwd=border.lwd, col="white")
    }
    axis(2, at=(1:(n))-0.5+y.offset, labels=rev(names), las=1, lty=0, hadj=0, line=name.line)
    if (show.range)
    {
        points((0:num.colors)/num.colors*n, rep(0, num.colors+1), col=colors, pch=19)
        axis(1, at=c(0,n/2, n), labels=as.character(signif(c(min(range), mean(range), max(range))),digits=3))
    }
}

make.diallel.skeleton <- function(num.strains, num.replicates)
{
    num.data <- num.replicates*num.strains^2

    # set up strain combinations
    mother.strain <- integer(num.data)
    father.strain <- integer(num.data)
    is.female <- integer(num.data)
    i <- 1
    for (m in 1:num.strains)
    {
        for (f in 1:num.strains)
        {
            for (r in 1:num.replicates)
            {
                mother.strain[i] <- m
                father.strain[i] <- f
                if (r/num.replicates > 0.5)
                {
                    is.female[i] <- 1
                }
                i <- i + 1
            }
        }
    }
    data <- data.frame(
            phenotype = rep(NA, num.data),
            father.strain.name=father.strain,
            mother.strain.name=mother.strain,
            is.female=is.female)
    data
}

expand.diallel.skeleton <- function(data,
        fill.missing=0,
        make.missing=0,
        phenotype="phenotype")
{
    required.cols <- c(phenotype, "is.female", "mother.strain.name", "father.strain.name")
    if (!all(required.cols %in% colnames(data)))
    {
        stop("Missing columns ", paste(collapse=", ", setdiff(required.cols, colnames(data))))
    }

    strain.name   <- sort(unique(c(as.character(data$mother.strain.name), as.character(data$father.strain.name))))
    data$mother.strain <- as.integer(factor(data$mother.strain.name, levels=strain.name))
    data$father.strain <- as.integer(factor(data$father.strain.name, levels=strain.name))

  data$original <- 1
    orig.mothers <- data$mother.strain
    orig.fathers <- data$father.strain
    orig.female  <- data$is.female
    if (0 < make.missing + fill.missing)
    {
        for (j in 1:length(strain.name))
        {
            for (k in 1:length(strain.name))
            {
                for (f in 0:1)
                {
                    num.wanted <- make.missing
                    if (0<fill.missing)
                    {
                        num.present <- sum(orig.mothers==j & orig.fathers==k & orig.female==f)
                        num.wanted <- max(0, fill.missing-num.present)+num.wanted
                    }
                    if (0==num.wanted) next
                    new.i <- (nrow(data)+1:num.wanted)
                    data[new.i,] <- NA
                    for (i in new.i)
                    {
                        data[i,c("mother.strain", "father.strain", "is.female", "original")] <- c(j,k,f,0)
                    }
                }
            }
        }
        data$mother.strain.name <- strain.name[data$mother.strain]
        data$father.strain.name <- strain.name[data$father.strain]
    }
    data$diallel.id <- 1:nrow(data)


    parentpair.name  <- paste(data$mother.strain.name, sep="_x_", data$father.strain.name)
    parentpair.id    <- as.integer(as.factor(parentpair.name))

    strainpair.name <- ifelse(data$mother.strain < data$father.strain,
                paste(data$mother.strain.name, sep="_x_", data$father.strain.name),
                paste(data$father.strain.name, sep="_x_", data$mother.strain.name))
    strainpair.id <- as.integer(as.factor(strainpair.name))

    parentpair.to.strainpair <- integer(unique(parentpair.id))
    parentpair.to.strainpair[parentpair.id] <- strainpair.id

    augmented.data <- cbind(data, data.frame(
            parentpair.name=parentpair.name,
            parentpair.id=parentpair.id,
            strainpair.name=strainpair.name,
            strainpair.id=strainpair.id,
            is.hybrid=ifelse(data$mother.strain==data$father.strain, 0, 1)))
    list(
            data=augmented.data,
            strain.name=strain.name,
            parentpair.to.strainpair=parentpair.to.strainpair,
            num.obs=nrow(data),
            num.strains=length(strain.name),
            num.parentpairs=length(parentpair.to.strainpair),
            num.strainpairs=length(unique(parentpair.to.strainpair))
            )
}

sample.until.converged <- function(jo, min.psrfq=1.1, max.rounds=10, burnin=1e3, ...)
{
    if (0<burnin)
    {
        update(jo, n.iter=burnin)
    }

    string.time <- function(a,b, divisor=1)
    {
        d <- difftime(b,a)
        paste(round(d/divisor,3), attr(d, "units"))
    }

    count     <- 0
    timer     <- list(Sys.time())
    mcmc      <- NULL
    converged <- FALSE
    while (!converged & count < max.rounds)
    {
        count <- count + 1
        cat("Sampling until converged: round ", count, "[", date(), "]:\n")
        mcmc <- coda.samples(jo, ...)
        timer[[count+1]] <- Sys.time()
        stat <- max(gelman.diag(mcmc)$psrf[,2])
        if (stat<min.psrfq)
        {
            cat("Converged to ", stat, "[")
            converged <- TRUE
        }
        else
        {
            cat("Unconverged to ", stat, "[")
        }
        cat("", string.time(timer[[count]],timer[[count+1]]), "]\n")
    }
    cat("Total time until convergence",
            string.time(timer[[1]], timer[[count+1]]),
            "[ average time of (", count, ") rounds =", string.time(timer[[1]], timer[[count+1]], count), "]\n")
    invisible(mcmc)
}

sim.diallel.effects <- function(data.list,
        var=list(a=0,d=0,m=0,
                strainpair=0, parentpair=0,
                delta.a=0, delta.d=0, delta.m=0,
                delta.strainpair=0, delta.parentpair=0,
                ind=1),
        beta=list(mu=0, hybrid=0, female=0, female.hybrid=0)
        )
{
    data <- data.list$data
    num.strains <- data.list$num.strains
    num.strainpairs <- data.list$num.strainpairs
    num.parentpairs <- data.list$num.parentpairs
    parentpair.to.strainpair <- data.list$parentpair.to.strainpair
    num.obs <- nrow(data.list$data)

    #----------------------------------------
    # specific combining ability effects
    # epistatic and ss-epistatic
    delta.strainpair <- rnorm(num.strainpairs, sd=sqrt(var$delta.strainpair))
    strainpair       <- rnorm(num.strainpairs, sd=sqrt(var$strainpair))

    # reciprocal epistatic and ss-reciprocal epistatic
    delta.parentpair <- rnorm(num.parentpairs,
            mean=delta.strainpair[parentpair.to.strainpair],
            sd=sqrt(var$delta.parentpair))
    parentpair <- rnorm(num.parentpairs,
            mean=strainpair[parentpair.to.strainpair],
            sd=sqrt(var$parentpair))

    #----------------------------------------
    # general combining ability effects

    delta.d <- rnorm(num.strains, mean=beta$female.hybrid, sd=sqrt(var$delta.d))
    d <- rnorm(num.strains, mean=beta$hybrid, sd=sqrt(var$d))

    delta.a <- rnorm(num.strains, sd=sqrt(var$delta.a))
    a <- rnorm(num.strains, sd=sqrt(var$a))

    delta.m <- rnorm(num.strains, sd=sqrt(var$delta.m))
    m <- rnorm(num.strains, sd=sqrt(var$m))

    ind <- rnorm(num.obs, sd=sqrt(var$ind))

    # record effects
    data.list$true.effects <- list(
        a=a, delta.a=delta.a,
        d=d, delta.d=delta.d,
        m=m, delta.m=delta.m,
        strainpair=strainpair, delta.strainpair=delta.strainpair,
        parentpair=parentpair, delta.parentpair=delta.parentpair)
    data.list$true.var     <- var
    data.list$true.beta    <- beta


    data <- data.list$data
    y.expect <- rep(NA, num.obs)
    for (i in 1:num.obs)
    {
        f <- ifelse(data$is.female[i], 0.5, -0.5)

        j <- data$mother.strain[i]
        aj.sex <- a[j] + f*delta.a[j]
        dj.sex <- d[j] + f*delta.d[j]
        mj.sex <- m[j] + f*delta.m[j]

        k <- data$father.strain[i]
        ak.sex <- a[k] + f*delta.a[k]
        dk.sex <- d[k] + f*delta.d[k]
        mk.sex <- m[k] + f*delta.m[k]

        mother.sex <- aj.sex + data$is.hybrid[i]*dj.sex + mj.sex
        father.sex <- ak.sex + data$is.hybrid[i]*dk.sex - mk.sex

        epi.jk.sex <- parentpair[data$parentpair.id[i]] + f*delta.parentpair[data$parentpair.id[i]]

        y.expect[i] <- beta$mu + mother.sex + father.sex + f*beta$female + epi.jk.sex
    }
    data$y.expect <- y.expect
    data$phenotype <- y.expect + ind

    data.list$data <- data
    data.list
}

timethis <- function(...)
{
    string.time <- function(a,b, divisor=1)
    {
        d <- difftime(b,a)
        paste(round(d/divisor,3), attr(d, "units"))
    }
    before <- Sys.time()
    out <- eval(...)
    after <- Sys.time()
    cat(string.time(before,after), "\n")
    out
}


draw.diallel.old <- function(data, phenotype="phenotype",
        main=paste("diallel", nrow(data), "obs"),
        show.text=FALSE,
        text.white.at=0.7,
        text.col=NA,
        text.fun=length
        )
{
    n <- max(data$mother.strain, data$father.strain)
    d.range <- range(data[,phenotype])

    strain.names <- 1:n
    if (!is.null(data$mother.strain.name))
    {
        nd <- unique(data[,c("mother.strain.name","mother.strain")])
        strain.names[nd$mother.strain] <- as.character(nd$mother.strain.name)
    }
    if (!is.null(data$father.strain.name))
    {
        nd <- unique(data[,c("father.strain.name","father.strain")])
        strain.names[nd$father.strain] <- as.character(nd$father.strain.name)
    }

    plot(c(0,n), c(0,n), axes=FALSE, type="n", ylab="", xlab="", main=main)
    title(ylab="father strain", line=1)
    title(xlab="mother strain", line=1)
    axis(1, at=(1:n)-0.5, labels=strain.names, lty=0, line=-1)
    if (!all(1:n==strain.names))
    {
        axis(1, at=(1:n)-0.5, labels=1:n, lty=0, line=-2)
    }
    axis(2, at=(1:n)-0.5, labels=rev(1:n), lty=0, las=1, line=-1)
    axis(3, at=(1:n)-0.5, labels=1:n, lty=0, line=-1)
    for (j in 1:n)
    {
        for (k in 1:n)
        {
            i <- data$mother.strain==k & data$father.strain==j
            d.cell <- data[i,phenotype]
            if (0==length(d.cell)) next
            d <- (mean(d.cell)-d.range[1])/diff(d.range)
            col <- gray(1-d)

            x.start <- j-1
            x.end   <- j
            y.start <- n-k
            y.end   <- n-k+1
            polygon(c(x.start,x.end,x.end,x.start), c(y.end, y.end, y.start, y.start), col=col, border=1, density=NA)

            if (show.text)
            {
                text.color <- text.col
                if (is.na(text.col))
                {
                    text.color <- gray(1-ifelse(d>text.white.at,0,1))
                }
                text(0.5*(x.start+x.end), 0.5*(y.start+y.end), labels=text.fun(d.cell), col=text.color)
            }
        }
    }
}

draw.diallel <- function(data, phenotype="phenotype",
        main=paste("diallel", nrow(data), "obs"),
        cex.strain.names=1,
        col=NULL,
        col.smoothness=2,
        data.matrix=NULL,
        data.range=range(na.rm=TRUE, data[,phenotype]),
        mothers.across=FALSE,
        na.action="na.exclude",
        na.symbol="cross",
        strain.names=NULL,
        strain.coords=NULL,
        show.strain.coords=TRUE,
        show.strain.names=TRUE,
        show.gridlines=TRUE,
        xlab=ifelse(mothers.across, "mother strain", "father strain"),
        ylab=ifelse(mothers.across, "father strain", "mother.strain")
        )
{
  if (!all(c("mother.strain", "father.strain") %in% colnames(data)))
  {
    msg="Need integer columns \"mother.strain\" and \"father.strain\" or factor columns \"mother.strain.name\" and \"father.strain.name\""
    if (!all(c("mother.strain.name", "father.strain.name") %in% colnames(data)))
    {
      stop(msg)
    }
    parents=c(as.character(data$mother.strain.name), as.character(data$father.strain.name))
    parents.int=as.integer(as.factor(parents))
    data$mother.strain=parents.int[1:nrow(data)]
    data$father.strain=parents.int[ (nrow(data)+1):(2*nrow(data)) ]
  }

    n <- max(data$mother.strain, data$father.strain)
    d.range <- range(data[,phenotype], na.rm=TRUE)
    if (!is.null(data.matrix))
    {
        stopifnot(identical(dim(data.matrix),c(n,n)))
    }
    else
    {
        data.matrix=matrix(nrow=n, ncol=n)
        for (j in 1:n)
        {
            for (k in 1:n)
            {
              #browser()
                i <- data$mother.strain==k & data$father.strain==j
                d.cell <- data[i,phenotype]
                if (0==length(d.cell))
                {
                  d <- NA
                }
                else if (all(is.na(d.cell)))
                {
                    d <- NA
                }
                else if (any(is.na(d.cell)) & "na.exclude"==na.action)
                {
                    d <- NA
                }
                else
                {
                    d <- mean(d.cell, na.rm=TRUE)
                }
                data.matrix[j,k] <- d
            }
        }
    }

    if (!mothers.across)
    {
      data.matrix=t(data.matrix)
    }

    if (is.null(col))
    {
        ncolors <- n^2*col.smoothness
        col=gray(1-(1:ncolors/ncolors))
    }

    image(0:n, 0:n, t(data.matrix[n:1,]), col=col, axes=FALSE, ylab="", xlab="", main=main,
        zlim=range(data.range, na.rm=TRUE))
    title(ylab=ylab, line=1.5)
    title(xlab=xlab, line=2)

    if (show.strain.coords)
    {
        if (is.null(strain.coords))
        {
            strain.coords <- 1:n
        }
        axis(1, at=(1:n)-0.5, labels=strain.coords, lty=0, line=-1)
        axis(2, at=(1:n)-0.5, labels=rev(strain.coords), lty=0, las=1, line=-0.5)
        axis(3, at=(1:n)-0.5, labels=strain.coords, lty=0, line=-0.5)
    }

    if (show.strain.names)
    {
        if (is.null(strain.names))
        {
            strain.names <- 1:n
            if (!is.null(data$mother.strain.name))
            {
                nd <- unique(data[,c("mother.strain.name","mother.strain")])
                strain.names[nd$mother.strain] <- as.character(nd$mother.strain.name)
            }
            if (!is.null(data$father.strain.name))
            {
                nd <- unique(data[,c("father.strain.name","father.strain")])
                strain.names[nd$father.strain] <- as.character(nd$father.strain.name)
            }
        }
        axis(1, at=(1:n)-0.5, labels=strain.names, lty=0, line=0, cex.axis=cex.strain.names)
    }

    if (show.gridlines) abline(h=0:n); abline(v=0:n)

    for (j in 1:n)
    {
        for (k in 1:n)
        {
            if (!is.na(data.matrix[j,k])) next
            if ("cross"==na.symbol)
            {
                segments(x0=rep(k-1,2), y0=(n-j)+c(0,1), x1=rep(k,2), y1=(n-j)+c(1,0))
            }
        }
    }
    invisible(data.matrix)
}


draw.strain.effects <- function(data,
        strain.names=rownames(data),
        strain.cols=rainbow(length(strain.names)),
        leftmargin=7, lwd=1)
{
    oldmar <- par("mar")
    on.exit(par(mar=oldmar))

    num.stats <- ncol(data)
    num.strains <- length(strain.names)
    par(mar=c(5,leftmargin,4,2))
    plot(range(data), c(1,num.stats), type="n", axes=FALSE, xlab="posterior mean effect", ylab="")
    axis(2, at=1:num.stats, labels=rev(colnames(data)), las=1)
    axis(1)
    strain.cols <- rainbow(num.strains)
    for (i in 1:num.strains)
    {
        lines(data[i,], num.stats:1, col=strain.cols[i], lwd=lwd)
    }
}



draw.diallel.sexed <- function(data, phenotype="phenotype", main=paste("diallel", nrow(data), "obs"))
{
    n <- max(data$mother.strain, data$father.strain)
    y.range <- range(data[,phenotype])

    strain.names <- 1:n
    if (!is.null(data$mother.strain.name))
    {
        nd <- unique(data[,c("mother.strain.name","mother.strain")])
        strain.names[nd$mother.strain] <- as.character(nd$mother.strain.name)
    }
    if (!is.null(data$father.strain.name))
    {
        nd <- unique(data[,c("father.strain.name","father.strain")])
        strain.names[nd$father.strain] <- as.character(nd$father.strain.name)
    }

    plot(c(0,n), c(0,2*n+1), axes=FALSE, type="n", ylab="", xlab="", main=main)
    title(ylab="father strain", line=1)
    title(xlab="mother strain", line=1)
    axis(1, at=(1:n)-0.5, labels=strain.names, lty=0, line=-1)
    axis(1, at=(1:n)-0.5, labels=1:n, lty=0, line=-2)
    axis(2, at=c((1:n)-0.5, (n+2):(2*n+1)-0.5), labels=rep(rev(1:n),2), lty=0, las=1, line=-1)
    text((1:n)-0.5, rep(n+0.5, n), labels=1:n)
    mtext("[females]", side=2, at=1.5*n, line=1)
    mtext("[males]", side=2, at=0.5*n, line=1)
    for (j in 1:n)
    {
        for (k in 1:n)
        {
            i <- data$mother.strain==k & data$father.strain==j
            for (is.female in 0:1)
            {
                si <- i & data$is.female==is.female
                y.cell <- data[si,phenotype]
                if (0==length(y.cell)) next
                x <- (mean(y.cell)-y.range[1])/diff(y.range)
                polygon(c(j-1,j,j,j-1), c(n-k+1, n-k+1, n-k, n-k)+is.female*(n+1), col=gray(1-x))
            }
        }
    }
}

prepare.jags.data <- function(data.list, model="sadme")
{
    components <- unlist(strsplit(model, ""))
    jags.data <- list(
            n = data.list$num.obs,
            phenotype = data.list$data$phenotype)

    if (any(components %in% c("s")))
    {
        jags.data$is.female <- as.integer(data.list$data$is.female)
    }
    if (any(components %in% c("d")))
    {
        jags.data$is.hybrid     <- data.list$data$is.hybrid
    }
    if (any(components %in% c("a", "d", "m")))
    {
        jags.data$num.strains   <- data.list$num.strains
        jags.data$mother.strain <- data.list$data$mother.strain
        jags.data$father.strain <- data.list$data$father.strain
    }
    if (any(components %in% c("e")))
    {
        jags.data$num.strainpairs           <- data.list$num.strainpairs
        jags.data$parentpair.to.strainpair  <- data.list$parentpair.to.strainpair

        jags.data$num.parentpairs           <- data.list$num.parentpairs
        jags.data$parentpair.id             <- data.list$parentpair.id
    }
    jags.data
}

read.sim.spec <- function(file)
{
    d <- read.table(file, header=FALSE)
    spec <- list(var=list(), beta=list())
    for (i in grep("^beta.*", d[,2]))
    {
        spec$beta[[sub("^beta.", "", d[i,2])]] <- as.numeric(d[i,1])
    }
    for (i in grep("^var.*", d[,2]))
    {
        spec$var[[sub("^var.", "", d[i,2])]] <- as.numeric(d[i,1])
    }
    spec
}

# GENERAL FUNCTIONS

ifow=function(test, yes, no){ if (test){return (yes)}; no }
