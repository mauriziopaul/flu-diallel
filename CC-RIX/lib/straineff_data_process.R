
get.column.index <- function(Y, col){
  which(colnames(Y) == col)
}

get.complete.row.index <- function(Y) {
  return(which(apply(Y, 1, function(x) sum(is.na(x))==0 )))
}

get.partial.row.index <- function(Y) {
  return(which(apply(Y, 1, function(x) sum(!is.na(x))>0 )))
}

## get phenotype data and other fixed effects
straineff.get.pgdata <- function(subjects,
                                 phenotype.file,
                                 phenotypecolname,
                                 subjectcolname,
                                 fixed.effects,
                                 random.effects){
  t=NULL
  sel=NULL
  X=NULL
  R=NULL
  ## second random effects
  R2=NULL

  Y <- read.csv(phenotype.file, header=T, row.names=NULL, sep="\t")

  suindex <- get.column.index(Y, subjectcolname)
  ## used for reading data in simulation study
  ## need to generate some fake phenotype

  if (is.null(phenotypecolname))
    phindex <- get.column.index(Y, subjectcolname)
  else{
    phindex <- unlist(strsplit(phenotypecolname,','))
    H <- length(phindex)
  }
  data <- Y[,phindex, drop=F]
  rownames(data) <- Y[, suindex]
  data <- data[subjects,,drop=F]
  sel <- names(get.partial.row.index(data))
  data <- data[sel, , drop=F]
  data <- t(data)
  N <- length(sel)

  X <- NULL
  if (!is.null(fixed.effects)){
    if (length(fixed.effects) == 1 && is.numeric(fixed.effects)){
      ##simulation
      num.dim<-fixed.effects
      X <-matrix(rbinom(N*num.dim, 1, 0.5), N, num.dim)
    }
    else {
      ## real data
      fixed <- Y[, fixed.effects, drop=F]
      if (!is.null(fixed.effects)){
        fixed <- as.matrix(fixed)
        rownames(fixed) <- Y[, suindex]
        fixed <- fixed[sel, , drop=F]
        for (i in 1:dim(fixed)[2]){
          if (is.numeric(fixed[, i])) {
            if (i==1){
              X <- fixed[, i]
            }
            else {
              X <- cbind(X,fixed[, i])
            }
          } else {
            x <- incidence.matrix(as.factor(fixed[, i]))
            if (i==1){
              X <- x[, -1]
            }
            else {
              X <- cbind(X, x[, -1])
            }
          }
        }
      }
      X <- as.matrix(X)
    }
  }

  if (!is.null(random.effects)){
    if (length(random.effects) <= 2) {
      random <- Y[, random.effects[1]]
      names(random) <- Y[, suindex]
      random <- random[subjects]
      random <- random[sel]
      R <- random
      if (length(random.effects) == 2) {
        random <- Y[, random.effects[2]]
        names(random) <- Y[, suindex]
        random <- random[subjects]
        random <- random[sel]
        R2 <- random
      }
    } else {
      stop("Cannot support multiple random effects now. ")
    }
  }
  list(phenotypes=data, sel=sel, X=X, R=R, R2=R2)
}



straineff.general.analysis <- function(model, data, subjects, phenotype.file,
                                       M=8, phenotypecolname, subjectcolname,
                                       transforms) {

  phenotypes <- pgdata$phenotypes
  sel <- pgdata$sel
  ## load fixed effects
  X <- pgdata$X

  ## load random effects
  K <- NULL
  Z <- NULL

  R <- pgdata$R
  R2 <- pgdata$R2

  ##transform phenotypes
  transforms <- unlist(strsplit(transforms,','))
  idx <- 1
  for (transform in transforms){
    if (transform=='identity') {
      phenotypes[idx,] <- phenotypes[idx,]
    }
    else if (transform=='log'){
      phenotypes[idx,] <- log10(phenotypes[idx,] + 0.01)
    }
    else if (transform=='sqrt'){
      phenotypes[idx,] <- sqrt(phenotypes[idx,])
    }
    else if (transform=='cubet'){
      # http://stackoverflow.com/questions/13236158/real-cube-root-of-a-negative-number-in-r
      phenotypes[idx,] <- sign(phenotypes[idx,]) * abs(phenotypes[idx,])^(1/3)
    }
    else if (transform=='sqr'){
      phenotypes[idx,] <- phenotypes[idx,]^2
    }
    else if (transform=='isqrt'){
      phenotypes[idx,] <- 1/sqrt(phenotypes[idx,])
    } else {
      stop("Unsupported tranformation function")
    }
    idx <- idx + 1
  }
  markers <- h$genotype$genome[which(h$genotype$genome$chromosome==chr), 1]
  h$strains <- h$full$strains
  data <- data[sel,]
  data[data < 0]=0
  data.deviation <- data[, (M + 1):(dim(data)[2])]
  data.additive <- met.or.vec(1,1)
  N <- dim(phenotypes)[2]
  data.additive[data.additive < 0] = 0
  data.additive <- data.additive[sel, ]
  ## TODO: if y is bionomial distribution
  method <- straineff.method.factory(model, data, data.additive[1:N, 1:M],
                                     data.deviation[1:N,], X,
                                     Y=as.matrix(phenotypes), Z, R, R2=R2)
  result <- method$estimate()
}
