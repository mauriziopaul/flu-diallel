run.diploffect.inla.through.genomecache <- function(formula, data, K=NULL,
                                                    genomecache, locus,
                                                    num.draws, use.dip.lincomb=TRUE, 
                                                    seed=1, gamma.rate=1, impute.on="SUBJECT.NAME", weights.on=NULL,  #"NUM.OBS", 
                                                    scale.response=TRUE,
                                                    do.augment.of.cache=FALSE){
  require(bagpipe.backend)
  #require(diplosoc.backend)
  
  
  weights.on <- weights.on[1]
  
  QTL.effect <- process.formula(formula, action="return.QTL.effect")
  if(QTL.effect == "additive") { add.only <- TRUE}
  if(QTL.effect == "full") { add.only <- FALSE }
  
  formula.string.new.data <- process.formula(formula, action="make.new.data.formula", impute.on=c(impute.on, weights.on))
  new.data <- model.frame(formula(formula.string.new.data), data=data)
  h <- bagpipe.backend::DiploprobReader$new(genomecache)
  founders <- h$getFounders()
  num.founders <- length(founders)
  # cache has rows per individual, not line
  if(!do.augment.of.cache){
    diplotype.matrix <- h$getLocusMatrix(locus=locus, model="full", subjects=new.data$SUBJECT.NAME)
    diplotype.matrix[diplotype.matrix < 0] <- 0
    # pick overlap between data and genomecache
    data <- data[as.character(data$SUBJECT.NAME) %in% h$getSubjects(),]
    new.data <- new.data[as.character(new.data$SUBJECT.NAME) %in% h$getSubjects(),]
  }
  # cache has rows per line
  if(do.augment.of.cache){
    diplotype.matrix <- h$getLocusMatrix(locus=locus, model="full", subjects=new.data[, impute.on])
    diplotype.matrix[diplotype.matrix < 0] <- 0
    
    rownames(diplotype.matrix) <- new.data$SUBJECT.NAME
    # pick overlap between data and genomecache
    data <- data[as.character(data[,impute.on]) %in% h$getSubjects(),]
    new.data <- new.data[as.character(new.data[,impute.on]) %in% h$getSubjects(),]
  }
  # pick overlap between data and genomecache
  diplotype.matrix <- diplotype.matrix[as.character(new.data$SUBJECT.NAME),]
  diploffect.inla.object <- run.diploffect.inla(formula=formula, data=data, K=K,
                                                num.founders=num.founders, prob.matrix=diplotype.matrix, add.only=add.only,
                                                num.draws=num.draws, use.dip.lincomb=use.dip.lincomb, 
                                                seed=seed, gamma.rate=gamma.rate,
                                                impute.on=impute.on, weights.on=weights.on, scale.response=scale.response, founders=founders, locus.name=locus)
  return(diploffect.inla.object)
}