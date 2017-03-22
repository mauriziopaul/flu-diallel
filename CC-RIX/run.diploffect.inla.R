run.diploffect.inla <- function(formula, data, K=NULL,
                                num.founders=8, prob.matrix, add.only,
                                num.draws, 
                                use.dip.lincomb=TRUE, seed=1, gamma.rate=1,
                                impute.on="SUBJECT.NAME", weights.on=NULL, #"NUM.OBS", 
                                scale.response=TRUE,
                                founders=NULL, locus.name=NULL){
  
  weights.on <- weights.on[1]
  
  formula.string.new.data <- process.formula(formula, action="make.new.data.formula", impute.on=c(impute.on, weights.on))
  fixed.formula.string <- process.formula(formula, action="make.fixed.formula")
  final.formula.string <- process.formula(formula, action="make.final.formula")
  new.data <- model.frame(formula(formula.string.new.data), data=data)
  rownames(new.data) <- new.data$SUBJECT.NAME

  # Setting up K structured random effect
  if(is.null(K)) { Z <- NULL }
  if(!is.null(K)){
    # taking overlap of data with kinship
    overlap <- intersect(new.data$SUBJECT.NAME, rownames(K))
    new.data <- new.data[overlap,]
    K <- K[overlap, overlap]
    
    Z <- factor(1:ncol(K))
  }
  
  Y <- new.data[,1]
  Y <- matrix(Y, 1, length(Y))
  rownames(Y) <- names(new.data)[1]
  colnames(Y) <- new.data$SUBJECT.NAME
  if(scale.response){
    Y[1,] <- scale(as.vector(Y), center=TRUE, scale=TRUE)
  }
  else{
    cat("Warning: Not centering and scaling the phenotype is not suggested.\n Doing so helps standardize hyperparameters across phenotypes\n")
  }
  X <- model.matrix(formula(fixed.formula.string), data=new.data)[,-1, drop=FALSE]
  rownames(X) <- new.data$SUBJECT.NAME
  if(dim(X)[2] == 0){
    X <- NULL
  }
  
  # Setting up heterscedasticity
  if(is.null(weights.on)){
    scale <- matrix(rep(1, nrow(Y)), 1, ncol(Y))
  }
  if(!is.null(weights.on)){
    scale <- matrix(new.data[, weights.on], 1, ncol(Y))
  }
  # Setting up sparse random effects
  Z2 <- NULL
  if(grepl(pattern="(1", x=final.formula.string, fixed=TRUE)){ # make random effect matrix to set up R2
    Z2.frame <- model.frame(formula(formula.string.new.data), data=data)
    random.var <- process.formula(formula, action="return.random.variable")
    Z2 <- Z2.frame[,random.var, drop=FALSE]
    rownames(Z2) <- Z2.frame$SUBJECT.NAME
  }
  
  # Founders
  if(is.null(founders)){
    founders <- LETTERS[1:num.founders]
  }
  
  # Prepping imputation in Diploffect
  imputation.map <- new.data[, c("SUBJECT.NAME", impute.on), drop=FALSE]
  imputation.map <- data.frame(SUBJECT.NAME=new.data$SUBJECT.NAME, impute.on=new.data[,impute.on])

  ROP <- NULL
  # if(do.ROP){
  #   ROP.ci <- fit.ROP(formula.string=final.formula.string, data=new.data, diplotype.matrix=dip.mat, K=K)
  #   
  #   lm.dip.fit <- fit.lm.ROP.dip(formula.string=final.formula.string, data=new.data, diplotype.matrix=dip.mat)
  #   lm.dip.ci <- extract.ROP.dip.ci(fit=lm.dip.fit, labels=colnames(dip.mat))
  #   lm.dos.fit <- fit.lm.ROP.dos(formula.string=formula.string, data=new.data, dosage.matrix=dip.mat.dos, diplotype.matrix=dip.mat)
  #   lm.ci <- extract.ROP.dos.ci(fit=lm.dos.fit, labels=c(colnames(dip.mat.dos), colnames(dip.mat)[-(1:M)]), M=M)
  #   
  #   lm.ci$diplotype.ci <- lm.dip.ci
  #   
  #   h2lmm.ci <- NULL
  #   if(!is.null(K)){
  #     h2lmm.dos.fit <- fit.h2lmm.ROP.dos(formula.string=final.formula.string, data=new.data, dosage.matrix=dip.mat.dos, diplotype.matrix=dip.mat, K=K)
  #     h2lmm.ci <- extract.ROP.dos.ci(fit=h2lmm.dos.fit, labels=c(colnames(dip.mat.dos), colnames(dip.mat)[-(1:M)]), M=M)
  #     h2lmm.dip.fit <- fit.h2lmm.ROP.dip(formula.string=formula.string, data=new.data, diplotype.matrix=dip.mat, K=K)
  #     h2lmm.dip.ci <- extract.ROP.dip.ci(fit=h2lmm.dip.fit, labels=colnames(dip.mat))
  #     
  #     h2lmm.ci$diplotype.ci <- h2lmm.dip.ci
  #   }
  #   
  #   ROP <- list(LM=lm.ci, H2LMM=h2lmm.ci)
  # }
  
  model <- ifelse(is.null(K), "inla", "inla.kinship")
  
  # Variance components
  random.var <- process.formula(formula, action="return.random.variable")
  if(length(random.var) > 0){ random.var <- paste0("random.", random.var) }
  var.components <- c(random.var, c("idx", "dom.idx", "poly.idx")[c(TRUE, !add.only, !is.null(K))])

  inla = INLAMethod$new()
  inla$init(model=model, data=prob.matrix, X=X, Y=Y, scale=scale, Z=Z, Z2=Z2, K=K)
  result = inla$estimate(num.draws=num.draws, 
                         family="gaussian",
                         num.threads=1,
                         use.dip.lincomb=use.dip.lincomb, 
                         gamma.rate=gamma.rate, 
                         this.seed=seed,
                         founder.names=founders,
                         ROP=ROP,
                         add.only=add.only,
                         var.components=var.components,
                         imputation.map=imputation.map)
  # Adding analysis information to Diploffect object
  genetic.effects <- c("additive", "dominant", "polygene")[c(TRUE, !add.only, !is.null(K))]
  result$analysis.id <- list(formula=final.formula.string, locus=locus.name, 
                             num.draws=num.draws, M=length(founders), founders=founders, 
                             genetic.effects=genetic.effects,
                             var.components=var.components)
  return(result)
}

process.formula <- function(formula, 
                            action=c("add.subjects", "return.effect", "make.final", "make.random.data.formula", "make.fixed.formula"),
                            impute.on="SUBJECT.NAME"){
  formula.string <- paste0(Reduce(paste, deparse(formula)))
  
  lh.formula <- trimws(unlist(strsplit(x=formula.string, split="~", fixed=TRUE)))[1]
  rh.formula <- trimws(unlist(strsplit(x=formula.string, split="~", fixed=TRUE)))[2]
  rh.formula <- trimws(unlist(strsplit(x=rh.formula, split="+", fixed=TRUE)))
  
  locus.index <- grepl(pattern="locus.", x=rh.formula, fixed=TRUE)
  random.effect.index <- grepl(pattern="\\(1\\s*\\|", x=rh.formula, perl=TRUE)
  
  if(action == "return.QTL.effect"){
    QTL.effect.type <- gsub(pattern="locus.", replacement="", x=rh.formula[locus.index])
    return(QTL.effect.type)
  }
  if(action == "make.fixed.formula"){
    fixed.formula.string <- paste(lh.formula, 
                                    paste(rh.formula[!(locus.index | random.effect.index)], collapse=" + "), 
                                    sep=" ~ ")
    return(fixed.formula.string)
  }
  if(action == "make.final.formula"){
    final.formula.string <- paste(lh.formula, 
                                  paste(rh.formula[!locus.index], collapse=" + "), 
                                  sep=" ~ ")
    return(final.formula.string)
  }
  if(action == "make.new.data.formula"){
    random.var <- trimws(gsub(pattern="\\(1\\s*\\|", replacement="", x=rh.formula[random.effect.index]))
    random.var <- gsub(pattern=")", replacement="", x=random.var)
    formula.string.new.data <- paste(lh.formula, 
                                  paste(paste(rh.formula[!(locus.index | random.effect.index)], collapse=" + "),
                                        paste(random.var, collapse=" + "), 
                                        paste(unique(c("SUBJECT.NAME", impute.on)), collapse=" + "), 
                                        sep=" + "),
                                  sep=" ~ ")
    return(formula.string.new.data)
  }
  if(action == "return.random.variable"){
    random.var <- trimws(gsub(pattern="\\(1\\s*\\|", replacement="", x=rh.formula[random.effect.index]))
    random.var <- gsub(pattern=")", replacement="", x=random.var)
    return(random.var)
  }
  if(action == "make.ROP.formula"){
    ROP.formula.string <- paste("scaled.y", 
                                paste(rh.formula, collapse=" + "), 
                                sep=" ~ ")
    return(ROP.formula.string)
  }
}

fit.ROP <- function(formula.string, data, diplotype.matrix, K){
  # Scaling and centering outcome
  data[,1] <- scale(data[,1], center=TRUE, scale=TRUE)
  names(data)[1] <- "scaled.y"
  ROP.formula.string <- process.formula(formula=formula(formula.string), action="make.ROP.formula")
  #M <- make.M(M=8)
  #data <- data.frame(data, diplotype.matrix %*% M)
  
  if(!is.null(K)){
    # h2lmm
    ROP.fit <- fit.h2lmm.ROP.dip(formula.string=ROP.formula.string, data=data, diplotype.)
  }
  else if(grepl(pattern="(1", x=ROP.formula.string, fixed=TRUE)){
    # lmer
    ROP.fit
  }
  else {
    # lm
    ROP.fit <- fit.lm.ROP.dip(formula.string=ROP.formula.string, data=data, diplotype.matrix=diplotype.matrix)
  }
  return(ROP.fit)
}

fit.lm.ROP.dip <- function(formula.string, data, diplotype.matrix){
  M <- make.M(M=36)
  rot.dip <- data.frame(diplotype.matrix %*% M)
  rotation.name <- paste("rotate.dip", 1:ncol(M), sep=".")
  names(rot.dip) <- rotation.name
  data <- data.frame(data, rot.dip)
  new.formula.string <- paste(formula.string,
                              paste(rotation.name, collapse=" + "),
                              sep=" + ")
  lm.fit <- lm(formula(new.formula.string), data=data)
  dip.effects <- M %*% coefficients(lm.fit)[grepl(pattern="rotate.dip", x=names(coefficients(lm.fit)), fixed=TRUE)]
  return(lm.fit)
}
fit.h2lmm.ROP.dip <- function(formula.string, data, diplotype.matrix, K){
  source("lmmbygls/scripts_to_source.R")
  h2lmm.fit <- lmmbygls(formula.string, data=data, K=K, use.par="h2")
  class(h2lmm.fit)[2] <- "lm"
  return(h2lmm.fit)
}
fit.lm.ROP.dos <- function(formula.string, data, dosage.matrix, diplotype.matrix){
  alt.formula.string <- paste(paste0("scaled.y ~ ", unlist(strsplit(formula.string, split="~"))[-1]), paste(colnames(dosage.matrix), collapse=" + "), paste(colnames(diplotype.matrix)[-(1:ncol(dosage.matrix))], collapse=" + "), sep=" + ")
  ROP.data <- data.frame(data, dosage.matrix, diplotype.matrix[,-(1:ncol(dosage.matrix))])
  ROP.data[,1] <- (ROP.data[,1] - mean(ROP.data[,1]))/sd(ROP.data[,1])
  names(ROP.data)[1] <- "scaled.y"
  lm.fit <- lm(alt.formula.string, data=ROP.data)
  return(lm.fit)
}
fit.h2lmm.ROP.dos <- function(formula.string, data, dosage.matrix, diplotype.matrix, K){
  source("lmmbygls/scripts_to_source.R")
  alt.formula.string <- paste(paste0("scaled.y ~ ", unlist(strsplit(formula.string, split="~"))[-1]), paste(colnames(dosage.matrix), collapse=" + "), paste(colnames(diplotype.matrix)[-(1:ncol(dosage.matrix))], collapse=" + "), sep=" + ")
  ROP.data <- data.frame(data, dosage.matrix, diplotype.matrix[,-(1:ncol(dosage.matrix))])
  ROP.data[,1] <- (ROP.data[,1] - mean(ROP.data[,1]))/sd(ROP.data[,1])
  names(ROP.data)[1] <- "scaled.y"
  h2lmm.fit <- lmmbygls(alt.formula.string, data=ROP.data, K=K, use.par="h2")
  class(h2lmm.fit)[2] <- "lm"
  return(h2lmm.fit)
}
extract.ROP.dos.ci <- function(fit, labels, M) {
  # dosages
  effects <- coef(fit)[labels[1:M]]
  nas <- is.na(effects)
  effects[nas] <- coef(fit)["(Intercept)"]
  effects[-nas] <- effects[!nas] + coef(fit)["(Intercept)"]

  ci.wide <- confint(fit, level=0.95)[labels[1:M],]
  ci.wide[nas,] <- confint(fit)["(Intercept)",]
  ci.wide[-nas,] <- ci.wide[!nas,] + confint(fit)["(Intercept)",]

  ci.narrow <- confint(fit, level=0.5)[labels[1:M],]
  ci.narrow[nas,] <- confint(fit)["(Intercept)",]
  ci.narrow[-nas,] <- ci.narrow[!nas,] + confint(fit)["(Intercept)",]
    
  colnames(ci.wide) <- c("quant.wide.left", "quant.wide.right")
  colnames(ci.narrow) <- c("quant.narrow.left", "quant.narrow.right")
  strain.ci <- list(med=effects, mu=effects, quant.narrow=ci.narrow, quant.wide=ci.wide)
  
  #deviations
  effects <- coef(fit)[labels[-(1:M)]]
  ci.wide <- confint(fit, level=0.95)[labels[-(1:M)],]
  ci.narrow <- confint(fit, level=0.5)[labels[-(1:M)],]
  
  colnames(ci.wide) <- c("quant.wide.left", "quant.wide.right")
  colnames(ci.narrow) <- c("quant.narrow.left", "quant.narrow.right")
  deviation.ci <- list(med=effects, mu=effects, quant.narrow=ci.narrow, quant.wide=ci.wide)
  return(list(strain.ci=strain.ci, deviation.ci=deviation.ci))
}
extract.ROP.dip.ci <- function(fit, labels) {
  effects <- coef(fit)[labels]
  nas <- is.na(effects)
  effects[nas] <- coef(fit)["(Intercept)"]
  effects[-nas] <- effects[!nas] + coef(fit)["(Intercept)"]
    
  ci.wide <- confint(fit, level=0.95)[labels,]
  ci.wide[nas,] <- confint(fit)["(Intercept)",]
  ci.wide[-nas,] <- ci.wide[!nas,] + confint(fit)["(Intercept)",]
    
  ci.narrow <- confint(fit, level=0.5)[labels,]
  ci.narrow[nas,] <- confint(fit)["(Intercept)",]
  ci.narrow[-nas,] <- ci.narrow[!nas,] + confint(fit)["(Intercept)",]
    
  colnames(ci.wide) <- c("quant.wide.left", "quant.wide.right")
  colnames(ci.narrow) <- c("quant.narrow.left", "quant.narrow.right")
  ci <- list(med=effects, mu=effects, quant.narrow=ci.narrow, quant.wide=ci.wide)
  return(ci)
}

run.diploffect.inla.summary.stats <- function(diploffect.object){
  deviation.ci <- deviation.mean <- diplotype.ci <- qtl.dom.ci <- SS.dom.ci <- SS.total.ci <- qtl.kinship.ci <- SS.poly.ci <- qtl.total.ci <- ROP <- nonqtl.ci.list <- nonqtl.ss.ci.list <- NULL
  M <- diploffect.object$analysis.id$M
  
  ##### Handling additional non-qtl variance components
  cat("Loading nonqtl proportion of variance explained (var components and sums of squares)... \n")
  var.components <- diploffect.object$analysis.id$var.components
  nonqtl.var.components <- var.components[grepl(var.components, pattern="random.")]
  if(length(nonqtl.var.components) > 0){
    nonqtl.ci.list <- nonqtl.ss.ci.list <- list()
    for(i in 1:length(nonqtl.var.components)){
      nonqtl.ci.list[[i]] <- load.nonqtl.h2.from.inla(diploffect.inla=diploffect.object, variable=nonqtl.var.components[i], scale=1)
      nonqtl.ss.ci.list[[i]] <- load.nonqtl.ss.h2.from.inla(diploffect.inla=diploffect.object, variable=nonqtl.var.components[i], scale=1)
    }
    names(nonqtl.ci.list) <- names(nonqtl.ss.ci.list) <- nonqtl.var.components
  }
  
  #### Additive
  if("additive" %in% diploffect.object$analysis.id$genetic.effects){
    strain.ci <- load.ci.from.inla(M=M, diploffect.inla=diploffect.object)
    cat("Loading strain effects...\n")
    qtl.add.ci <- load.qtl.h2.from.inla(diploffect.inla=diploffect.object, effect="add", scale=1)
    cat("Loading additive heritability...\n")
    SS.add.ci <- load.SS.h2.from.inla(diploffect.inla=diploffect.object, effect="add", scale=1)
    cat("Loading additive sums of squares...\n")
  }
  
  ##### Dominant
  if("dominant" %in% diploffect.object$analysis.id$genetic.effects){
    deviation.ci <- load.deviated.ci.from.inla(M, diploffect.inla=diploffect.object, scale=1)
    cat("Loading deviation effects...\n")
    deviation.mean <- load.deviated.mean.from.inla(M, diploffect.inla=diploffect.object, scale=1)
    cat("Loading deviation means...\n")
    diplotype.ci <- load.diplotypes.ci.from.inla(M, diploffect.inla=diploffect.object, scale=1)
    cat("Loading diplotype means...\n")
    qtl.dom.ci <- load.qtl.h2.from.inla(diploffect.inla=diploffect.object, effect="dom", scale=1)
    cat("Loading dominant heritability...\n")
    qtl.total.ci <- load.total.qtl.h2.from.inla(diploffect.inla=diploffect.object, scale=1)
    cat("Loading total QTL heritability...\n")
    SS.dom.ci <- load.SS.h2.from.inla(diploffect.inla=diploffect.object, effect="dom", scale=1)
    cat("Loading dominant sums of squares...\n")
    SS.total.ci <- load.SS.h2.from.inla(diploffect.inla=diploffect.object, effect="qtl", scale=1)
    cat("Loading total QTL sums of squares...\n")
  }
  if("polygene" %in% diploffect.object$analysis.id$genetic.effects){
    qtl.kinship.ci <- load.qtl.h2.from.inla(diploffect.inla=diploffect.object, effect="kinship", scale=1)
    cat("Loading polygene heritability...\n")
    SS.poly.ci <- load.SS.h2.from.inla(diploffect.inla=diploffect.object, effect="poly", scale=1)
    cat("Loading polygene sums of squares...\n")
  }
  if(!is.null(diploffect.object$ROP)){
    ROP <- load.ROP.estimates(diploffect.inla=diploffect.object)
  }
  
  ci.results <- list(strain.ci=strain.ci, deviation.ci=deviation.ci, deviation.mean=deviation.mean, diplotype.ci=diplotype.ci,
                     qtl.total.ci=qtl.total.ci, qtl.add.ci=qtl.add.ci, qtl.dom.ci=qtl.dom.ci, kinship.ci=qtl.kinship.ci, 
                     SS.total.ci=SS.total.ci, SS.add.ci=SS.add.ci, SS.dom.ci=SS.dom.ci, SS.poly.ci=SS.poly.ci,
                     ROP=ROP, nonqtl.ci.list=nonqtl.ci.list, nonqtl.ss.ci.list=nonqtl.ss.ci.list,
                     analysis.id=diploffect.object$analysis.id)
  return(ci.results)
}


make.M <- function(M){
  j <- M
  k <- (-1 + sqrt(j))*(j - 1)^(-3/2)
  m <- 1/sqrt(j - 1)
  M <- diag(j - 1)
  M[M == 0] <- -k
  M <- rbind(M, rep(-m, j - 1))
  return(M)
}