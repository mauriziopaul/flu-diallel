#library('glmnet')

StrainEffMethod <- setRefClass("StrainEffMethod",
  fields = list(X_="matrix", Y_="matrix", scale_="matrix", model_="character",
    data_="matrix", N_="numeric",
    Z_="matrix", 
    Z2_="list",
    K_="matrix", cholK_="matrix", rawZ_="factor", rawZ2_="data.frame"),
  methods = list(
    init = function(model, data, X, Y, scale, Z, Z2, K) {
      'TODO: call by data.additive[1:N,1:8]'
      Y_ <<- Y
      if (!is.null(X)) X_ <<- X
      data_ <<- data
      model_ <<- model

      N_ <<- dim(Y)[2]

      if (!is.null(Z)) {
        Z_ <<- incidence.matrix(as.factor(Z))
        rawZ_ <<- Z
      }
      if (!is.null(Z2)) {
        Z2_ <<- list()
        for(i in 1:ncol(Z2)){
          Z2_[[i]] <<- incidence.matrix(as.factor(Z2[,i]))
        }
        names(Z2_) <<- names(Z2)
        rawZ2_ <<- Z2
      }
      if (!is.null(K)) {
        K_ <<- K
      }
    },
    can.estimate = function(model) {
      'check whether the class can be used for estimation'
      stop("Not Implemented")
    },
    estimate = function() {
      'Estimate strain effects'
      stop("Not Implemented")
    },
    get.straineff = function(M, result) {
      'Estimate strain effects'
      stop("Not Implemented")
    },
    is.straineff = function() {
      'Is strain effects model?'
      return(FALSE)
    },
    is.straineff.mcmc = function() {
      'Is strain effects model with MCMC sampling?'
      return(FALSE)
    },
    has.deviation.effects = function() {
      'Does the method supply deviation estimation?'
      return (FALSE)
    },
    get.deviation.effects = function(M, result) {
      'Estimate strain effects'
      stop("Not Implemented")
    },
    has.deviation.effects = function() {
      return(FALSE)
    },
    get.diplotype.effects = function(M, result) {
      beta = get.straineff(M, result)
      if (!has.deviation.effects()) {
        deviation.effects <- mat.or.vec(1, M * (M + 1) / 2)
      }
      else {
        deviation.effects <- get.deviation.effects(M, result)
      }
      return(calculate.diplotype.effects(beta, deviation.effects))
    },
    get.haplotype.score = function(result) {
      return (0)
    },
    has.haplotype.score = function(result) {
      return(F)
    },
    get.display.name = function(result) {
      return(model_)
    },
    get.extra = function(result) {
      return (NULL)
    }
  )
)

get.straineff.method <- function(model) {
  #candidate.methods = c(GLMnetMethod, LMMethod, SurrogateMethod, StrainEffModel, INLAMethod)
  candidate.methods = c(INLAMethod)
  for (m in candidate.methods) {
    ret = m$new()
    if (ret$can.estimate(model)) {
      return (ret)
    }
  }
  return (NULL)
}

straineff.method.factory <- function(model, data, data.additive, data.deviated,
                                     X, Y, scale, R, R2, K=NULL) {
  method = get.straineff.method(model)
  if (!is.null(method)) {
    method$init(model, data, data.additive, data.deviated, X, Y, scale, R, R2, K)
    return (method)
  }
  else {
    stop("Did you forget to put your method file in straineff.method.R?")
  }
}
