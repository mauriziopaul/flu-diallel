source("lib/mcmc_helper.R")
source("lib/library_diallel_helper.R")

get.ci.from.mcmc <- function(dat, wanted) {
  dat[,wanted] = dat[,wanted]- rowMeans(dat[,wanted])
  which.wanted=ifow(is.integer(wanted), wanted, match(wanted, varnames(dat)))
  num.wanted=length(which.wanted)
  chain <- mcmc.stack(dat)
  mu    <- colMeans(chain[,which.wanted])
  med   <- apply(coda::HPDinterval(chain, prob=0.01)[which.wanted,],
                 1, mean)
  prob.wide=0.95
  prob.narrow=0.50
  hpd.wide    <- coda::HPDinterval(chain, prob=prob.wide)[which.wanted,]
  hpd.narrow  <- coda::HPDinterval(chain, prob=prob.narrow)[which.wanted,]
  return (list(med=med, mu=mu, hpd.narrow=hpd.narrow, hpd.wide=hpd.wide))
}

load.ci.from.mcmc <- function(M, file, H=1, h=1) {
  load(file)
  mcmc = result[[1]][[1]]
  extra <- straineff.extra(H, h)
  ci <- get.ci.from.mcmc(mcmc, wanted=paste('beta[', extra, 1:M, ']', sep=""))
  ci
}

plot.ci <- function(midvals, narrow.intervals, wide.intervals,
                    names=1:length(midvals),
                    add=FALSE,
                    xlab="Estimate",
                    xlab.line=2.5,
                    ylab="",
                    yaxis=TRUE,
                    name.margin=6,
                    name.line=4,
                    pch.midvals=19,
                    col=rep("black", length(midvals)),
                    col.midvals=col,
                    shift=0,
                    type="p",
                    use.this.lim=NULL,
                    main="",
                    ...)
{
  nvals <- length(midvals)
  col.midvals <- rep(col.midvals, length.out=nvals)
  y.pos <- (1:nvals)-0.5
  if(!add){
    if(is.null(use.this.lim)){
      lim <- range(c(wide.intervals,narrow.intervals,midvals), na.rm=TRUE)
      lim <- c(-1,1) * diff(lim)*0.1 + lim 
    }
    if(!is.null(use.this.lim)){
      lim <- use.this.lim
    }
    
    mar <- c(5, name.margin, 4, 2)+0.1
    oldmar <- par(mar=mar); on.exit(par(mar=oldmar))
    plot(lim, c(0,nvals+0.5), type="n", axes=FALSE, ylab=ylab, xlab="", main=main, ...)
    title(xlab=xlab, line=xlab.line)
    axis(1)
    axis(3, line=-1)
    if(yaxis){
      axis(2, at=y.pos, labels=rev(names), las=1, lty=0, hadj=0, line=name.line)
    }
  }
  if("p"==type){
    for(i in 1:nvals){
      pos <- nvals-i + 0.5 + shift
      lines(wide.intervals[i,], rep(pos,2), col=col[i])
      lines(narrow.intervals[i,], rep(pos,2), lwd=3, col=col[i])
      points(midvals[i], pos, pch=pch.midvals, col=col.midvals[i])
    }
  }
  invisible(rev(y.pos))
}

prepare.additive.dominant.ratio.posterior <- function(files, num.draws=1000) {
  samples = c()
  weights = c()
  if (!is.vector(files)) {
    files = c(files)
  }
  for (file in files) {
    load(file)
    print(file)
    logmliks = as.vector(scale(result$results$logmliks, scale=F))
    mliks = exp(logmliks)
    ratio = lapply(1:length(mliks), function (i) {
      hyper.samples = result$results$hyper.samples[[i]]
      ## bug in inla, this is not log scale
      precision.add = hyper.samples[, 'Log precision for idx1 in user-scale']
      precision.dom = hyper.samples[, 'Log precision for dom.idx in user-scale']

      (1 / (precision.add ^ 2)) / (1 / (precision.add ^ 2) + 1 / (precision.dom ^2))
    })

    samples = as.vector(unlist(ratio))
    weights = rep(mliks, each=num.draws)
  }
  weights = weights / sum(weights)
  density(samples, weights=weights, from=0, to=1, kernel='cosine')
}

draw.diallel <- function(data.matrix, sn=hs.sn, main="") {
  data.range=range(as.vector(data.matrix))
  n = dim(data.matrix)[1]
  col.smoothness = 2
  ncolors <- n^2*col.smoothness
  col=gray(1-(1:ncolors/ncolors))
  image(0:n, 0:n, t(data.matrix[n:1,]), col=col, axes=FALSE, ylab="", xlab="", main=main,
        zlim=range(data.range, na.rm=TRUE))

  strain.coords <- 1:n
  axis(1, at=(1:n)-0.5, labels=strain.coords, lty=0, line=-0.5)
  axis(2, at=(1:n)-0.5, labels=rev(strain.coords), lty=0, las=1, line=-0.5)
  axis(3, at=(1:n)-0.5, labels=strain.coords, lty=0, line=-0.5)

  axis(1, at=(1:n)-0.5, labels=sn, lty=0, line=0.5, cex.axis=0.8)
  abline(h=0:n); abline(v=0:n)
}

plot.straineff.ci <- function(inla.diploffect.ci, sn=NULL, xlab="Haplotype Effects", main="", flip=TRUE, ...) {
  ci <- inla.diploffect.ci$strain.ci
  if(is.null(sn)){
    sn <- inla.diploffect.ci$analysis.id$founders
  }
  
  if(flip){ order <- rev(1:nrow(ci$quant.narrow)) }
  if(!flip){ order <- 1:nrow(ci$quant.narrow) }
  main <- c(paste(inla.diploffect.ci$analysis.id$formula, paste0("locus(", inla.diploffect.ci$analysis.id$locus, ")"), sep=" + "),
            paste("INLA samples:", inla.diploffect.ci$analysis.id$num.draws))
  ypos <- plot.ci(ci$med[order], ci$quant.narrow[order,], ci$quant.wide[order,], names=sn[order],
                  xlab=xlab, col.midvals="white",
                  pch.midvals="|", type="p", main=main, ...)
  points(ci$mu[order], ypos, pch="|")
  abline(v=0, lty=2)
}

plot.deviation.ci <- function(inla.diploffect.ci, sn=NULL, xlab="Dominant Deviation Effects", flip=TRUE, ...) {
  ci <- inla.diploffect.ci$deviation.ci
  if(is.null(sn)){
    founders <- inla.diploffect.ci$analysis.id$founders
    full.to.dosages <- t(straineff.mapping.matrix.happy(M=length(founders)))
    
    sn <- apply(full.to.dosages[-(1:8),], 1, function(x) paste(founders[sort(which(x==1), decreasing=FALSE)], collapse=" x "))
  }
  if(flip){ order <- rev(1:nrow(ci$quant.narrow)) }
  if(!flip){ order <- 1:nrow(ci$quant.narrow) }
  main <- c(paste(inla.diploffect.ci$analysis.id$formula, paste0("locus(", inla.diploffect.ci$analysis.id$locus, ")"), sep=" + "),
            paste("INLA samples:", inla.diploffect.ci$analysis.id$num.draws))
  ypos <- plot.ci(ci$med[order], ci$quant.narrow[order,], ci$quant.wide[order,], names=sn[order],
                  xlab=xlab, col.midvals="white",
                  pch.midvals="|", type="p", main=main, ...)
  points(ci$mu[order], ypos, pch="|")
  abline(v=0, lty=2)
}

plot.diplotype.ci <- function(inla.diploffect.ci, sn=NULL, xlab="Diplotype Effects", flip=TRUE, ...) {
  ci <- inla.diploffect.ci$diplotype.ci
  if(is.null(sn)){
    founders <- inla.diploffect.ci$analysis.id$founders
    full.to.dosages <- t(straineff.mapping.matrix.happy(M=length(founders)))
    
    sn <- c(paste(founders, founders, sep=" x "),
            apply(full.to.dosages[-(1:8),], 1, function(x) paste(founders[sort(which(x==1), decreasing=FALSE)], collapse=" x ")))
  }
  if(flip){ order <- rev(1:nrow(ci$quant.narrow)) }
  if(!flip){ order <- 1:nrow(ci$quant.narrow) }
  main <- c(paste(inla.diploffect.ci$analysis.id$formula, paste0("locus(", inla.diploffect.ci$analysis.id$locus, ")"), sep=" + "),
            paste("INLA samples:", inla.diploffect.ci$analysis.id$num.draws))
  par(mar=c(3.2,3.2,2,1), mgp=c(2.2,.7,0), tck=-.01, las=1)
  ypos <- plot.ci(ci$med[order], ci$quant.narrow[order,], ci$quant.wide[order,], names=sn[order],
                  xlab=xlab, col.midvals="white",
                  pch.midvals="|", type="p", main=main, ...)
  points(ci$mu[order], ypos, pch="|")
  abline(v=0, lty=2)
}

plot.comparison.cis <- function(ci.list, analysis.id, labels=NULL, sn, xlab="Haplotype Effects", add.numbers=FALSE,
                                use.this.lim=NULL, comp.col=c("black", "red"), ...) 
  {
  par(mar=c(3.2,3.2,2,1), mgp=c(2.2,.7,0), tck=-.01, las=1)
  main <- c(paste(analysis.id$formula, paste0("locus(", analysis.id$locus, ")"), sep=" + "),
            paste("INLA samples:", analysis.id$num.draws))
  add = FALSE
  step = 0.4 / (length(ci.list) - 1)
  shift = 0.2 ## start
  if(is.null(use.this.lim)){
    lim <- range(ci.list, na.rm=TRUE)
    lim <- c(-1,1) * diff(lim)*0.1 + lim 
  }
  else{
    lim <- use.this.lim
  }
  for(i in 1:length(ci.list)){
    ci <- ci.list[[i]]
    if(add){
      ypos <- plot.ci(ci$med, ci[[3]], ci[[4]], names=sn, add=add,
                      xlab=xlab, col.midvals="white", col=rep(comp.col[i], length(ci$med)),
                      pch.midvals="|", type="p", shift=shift, main=main, use.this.lim=lim, ...)
    }
    if(!add){
      ypos <- plot.ci(ci$med, ci[[3]], ci[[4]], names=sn, add=add,
                      xlab=xlab, col.midvals="white", col=rep(comp.col[i], length(ci$med)),
                      pch.midvals="|", type="p", shift=shift, use.this.lim=lim, main=main, ...)
    }
    points(ci$mu, ypos + shift, pch="|")
    if(add.numbers){
      mtext(text=paste(round(ci$mu, 2), paste0("(", round(ci[[4]][, 1], 2), ", ", round(ci[[4]][, 2], 2), ")")), side=4, at=ypos + shift, adj=1, col=comp.col[i])
    }
    if(!add){ add <- TRUE }
    shift <- shift - step
  }
  abline(v=0, lty=2)
  
  if(!is.null(labels)){
    legend('bottomleft', labels,
           col=1:length(labels), lty=1, bty='n', lwd=3)
  }
}

plot.qtl.varexp.ci <- function(inla.diploffect.ci, xlab="QTL Effects (Variance Explained)", ...){
  join.ci <- function(ci.list) {
    combine.ci <- list()
    med <- NULL
    mu <- NULL
    quant.narrow <- NULL
    quant.wide <- NULL
    for (i in 1:length(ci.list)) {
      med <- c(med, ci.list[[i]]$med)
      mu <- c(mu, ci.list[[i]]$mu)
      quant.narrow <- rbind(quant.narrow, ci.list[[i]]$quant.narrow)
      quant.wide <- rbind(quant.wide, ci.list[[i]]$quant.wide)
    }
    combine.ci <- list(med=med, mu=mu, quant.narrow=quant.narrow, quant.wide=quant.wide)
    return(combine.ci)
  }
  
  SS.ci.list <- list(inla.diploffect.ci$SS.add.ci)
  h2.ci.list <- list(inla.diploffect.ci$qtl.add.ci)
  
  effect.labels <- c("QTL Additive")
  if("dominant" %in% inla.diploffect.ci$analysis.id$genetic.effects){
    SS.ci.list[[length(SS.ci.list)+1]] <- inla.diploffect.ci$SS.dom.ci
    SS.ci.list[[length(SS.ci.list)+1]] <- inla.diploffect.ci$SS.total.ci
    
    h2.ci.list[[length(h2.ci.list)+1]] <- inla.diploffect.ci$qtl.dom.ci
    h2.ci.list[[length(h2.ci.list)+1]] <- inla.diploffect.ci$qtl.total.ci
    
    effect.labels <- c(effect.labels, "QTL Dominant", "QTL Total")
  }
  if("polygene" %in% inla.diploffect.ci$analysis.id$genetic.effects){
    SS.ci.list[[length(SS.ci.list)+1]] <- inla.diploffect.ci$SS.poly.ci
    h2.ci.list[[length(h2.ci.list)+1]] <- inla.diploffect.ci$kinship.ci
    
    effect.labels <- c(effect.labels, "QTL Polygene")
  }
  
  SS.ci <- join.ci(SS.ci.list)
  h2.ci <- join.ci(h2.ci.list)
  
  plot.comparison.cis(ci.list=list(SS.ci, h2.ci), 
                      analysis.id=inla.diploffect.ci$analysis.id, 
                      labels=c("SS", "VC"), sn=effect.labels, xlab=xlab)
  
}

plot.varexp.ci <- function(inla.diploffect.ci, xlab="Variance Explained", add.numbers=FALSE, ...){
  join.ci <- function(ci.list) {
    combine.ci <- list()
    med <- NULL
    mu <- NULL
    quant.narrow <- NULL
    quant.wide <- NULL
    for (i in 1:length(ci.list)) {
      med <- c(med, ci.list[[i]]$med)
      mu <- c(mu, ci.list[[i]]$mu)
      quant.narrow <- rbind(quant.narrow, ci.list[[i]]$quant.narrow)
      quant.wide <- rbind(quant.wide, ci.list[[i]]$quant.wide)
    }
    combine.ci <- list(med=med, mu=mu, quant.narrow=quant.narrow, quant.wide=quant.wide)
    return(combine.ci)
  }
  
  SS.ci.list <- list()
  h2.ci.list <- list()
  
  effect.labels <- NULL
  if(!is.null(inla.diploffect.ci$nonqtl.ci.list)){
    for(i in 1:length(inla.diploffect.ci$nonqtl.ci.list)){
      h2.ci.list[[i]] <- inla.diploffect.ci$nonqtl.ci.list[[i]]
      SS.ci.list[[i]] <- inla.diploffect.ci$nonqtl.ss.ci.list[[i]]
      effect.labels <- c(effect.labels, names(inla.diploffect.ci$nonqtl.ci.list)[i])
    }
  }

  SS.ci.list[[length(SS.ci.list)+1]] <- inla.diploffect.ci$SS.add.ci
  h2.ci.list[[length(h2.ci.list)+1]] <- inla.diploffect.ci$qtl.add.ci
  effect.labels <- c(effect.labels, "QTL Additive")
  if("dominant" %in% inla.diploffect.ci$analysis.id$genetic.effects){
    SS.ci.list[[length(SS.ci.list)+1]] <- inla.diploffect.ci$SS.dom.ci
    SS.ci.list[[length(SS.ci.list)+1]] <- inla.diploffect.ci$SS.total.ci
    
    h2.ci.list[[length(h2.ci.list)+1]] <- inla.diploffect.ci$qtl.dom.ci
    h2.ci.list[[length(h2.ci.list)+1]] <- inla.diploffect.ci$qtl.total.ci
    
    effect.labels <- c(effect.labels, "QTL Dominant", "QTL Total")
  }
  if("polygene" %in% inla.diploffect.ci$analysis.id$genetic.effects){
    SS.ci.list[[length(SS.ci.list)+1]] <- inla.diploffect.ci$SS.poly.ci
    h2.ci.list[[length(h2.ci.list)+1]] <- inla.diploffect.ci$kinship.ci
    
    effect.labels <- c(effect.labels, "QTL Polygene")
  }
  
  SS.ci <- join.ci(SS.ci.list)
  h2.ci <- join.ci(h2.ci.list)
  
  plot.comparison.cis(ci.list=list(SS.ci, h2.ci), 
                      analysis.id=inla.diploffect.ci$analysis.id, 
                      labels=c("SS", "VC"), sn=effect.labels, xlab=xlab, add.numbers=add.numbers, use.this.lim=c(0, 1))
  abline(v=1, lty=2)
  
}

plot.diallel <- function(inla.diploffect.ci, sn=NULL) {
  if(is.null(sn)){
    sn <- inla.diploffect.ci$analysis.id$founders
  }
  diplotype <- inla.diploffect.ci$diplotype.ci$mu
  M <- inla.diploffect.ci$analysis.id$M
  main <- c(paste(inla.diploffect.ci$analysis.id$formula, paste0("locus(", inla.diploffect.ci$analysis.id$locus, ")"), sep=" + "),
            paste("INLA samples:", inla.diploffect.ci$analysis.id$num.draws))
  par(mar=c(3.5, 3.2, 4.5, 1.0))
  mapping = straineff.mapping.matrix(M, 'happy')
  data.matrix = matrix(0, M, M)
  for (i in 1:ncol(mapping)) {
    matrix.idx <- which(mapping[,i] != 0)
    if (length(matrix.idx) == 1) {
      data.matrix[matrix.idx, matrix.idx] <- diplotype[i]
    }
    else if (length(matrix.idx) == 2) {
      data.matrix[matrix.idx[1], matrix.idx[2]] <- data.matrix[matrix.idx[2], matrix.idx[1]] <- diplotype[i]
    }
  }
  draw.diallel(data.matrix, sn=sn, main=main)
}
