require(INLA)
get.combined.marginal = function(samples, fieldname, rowidx, mliks) {
  S = length(samples)
  stopifnot(S == length(mliks))
  s = rep(0, S)
  get.max.from.marginal <- function(s) {
    return(max(s$marginals.random[[fieldname]][[rowidx]][, 1]))
  }
  get.min.from.marginal <- function(s) {
    return(min(s$marginals.random[[fieldname]][[rowidx]][, 1]))
  }
  left = min(unlist(lapply(samples, get.max.from.marginal)))
  right = max(unlist(lapply(samples, get.min.from.marginal)))

  x = seq(left, right, length.out=200)
  get.weighted.density <- function(x) {
    s = rep(0, S)
    for (i in 1:S) {
      s[i] = inla.dmarginal(
         x, samples[[i]]$marginals.random[[fieldname]][[rowidx]])
    }
    return (weighted.mean(s, mliks))
  }
  y = unlist(lapply(x, get.weighted.density))
  return (data.frame(x, y))
}

span.matrix <- function(x, pattern, H) {
  stopifnot(dim(x)[2] == sum(pattern))
  nx = matrix(0, dim(x)[1] * H, dim(x)[2] * H)
  old.col.start = 1
  col.start = 1
  for (p in pattern) {
    row.start = 1
    for (h in 1:H) {
      print(nx[row.start:(row.start + dim(x)[1] - 1), col.start:(col.start + p - 1)])
      nx[row.start:(row.start + dim(x)[1] - 1), col.start:(col.start + p - 1)] = x[, old.col.start:(old.col.start + p - 1)]
      row.start = row.start + dim(x)[1]
      col.start = col.start + p
    }
    old.col.start = old.col.start + p
  }
  return (nx)
}

INLAMethod <- setRefClass("INLAMethod",
  contains = "StrainEffMethod",
  methods = list(
    can.estimate = function(model) {
      model_ <<- model
      return (grepl('inla', model))
    },
    use.weight = function() {
      return (!grepl('noweight', model_))
    },
    use.kinship = function() {
      return (grepl('kinship', model_))
    },
    estimate = function(num.draws=200, family="gaussian", num.threads=4) {
      S = num.draws
      ## The number of phenotypes
      H = dim(Y_)[1]
      results = list(logmliks=rep(0, S), betas = matrix(0, S, M_ * H),
        summary.random=list(), inla.objects=list(), hyper.samples=list())
      T = M_ * (M_ + 1) / 2
      diplotypes = list()
      case = 1
      while(TRUE) {

        s = dim(data_)[2]
        data = matrix(0, N_, s)

        diplotype = rep(0, N_)
        mapping = straineff.mapping.matrix(M_, 'happy')
        ## Impute the diplotype (haplotype pair) status
        x = matrix(0, N_, M_)
        for (i in 1:N_) {
          data[i, ] = rmultinom(1, 1, data_[i, ])
          idx = which(data[i, ] == 1)
          x[i, ] = mapping[, idx]
          diplotype[i] = idx
        }
        diplotypes[[case]] = diplotype
        ## If the diplotype does not exist in the imputed data,
        ## remove the corresponding diplotype effects.
        selected = c()
        for (i in (M_ + 1):T) {
          if (sum(data[, i]) != 0) {
            selected = c(selected, i)
          }
        }
        dom.M <- length(selected)

        ## the following script is based on the example at the bottom
        ## of http://www.r-inla.org/models/tools
        if (dim(X_)[1]==0) {
          x = cbind(rep(1, N_), x)
        }
        else {
          x = cbind(rep(1, N_), X_, x)
        }
        ##print(qr(x)$rank)
        if (dom.M != 0){
          x = cbind(x, data[,selected])
        }

        random.M <- 0
        if (dim(R_)[1] != 0) {
          if (use.kinship()) {
            print("warning, Kinship function in INLA is too slow")
            ck = K_
            x = cbind(x, ck)
            random.M = dim(ck)[1]
          }
          else {
            x = cbind(x, R_)
            random.M = dim(R_)[2]
          }
        }

        random2.M <- 0
        if (dim(R2_)[1] != 0) {
          x = cbind(x, R2_)
          random2.M = dim(R2_)[2]
        }

        fixed.M <- dim(X_)[2]
        fixed =NA
        ## kronecker product
        if (H == 2) {
          x = diag(H) %x% x
        }
        if (H == 1) {
          y = Y_[1, ]
          data <- list(y=y,
                       intercept=c(1, rep(NA, fixed.M + M + dom.M + random.M + random2.M)),
                       fixed=c(NA, 1:fixed.M, rep(NA,  M + dom.M + random.M + random2.M)),
                       idx1=c(rep(NA, 1 +  fixed.M), 1:M, rep(NA, dom.M + random.M + random2.M)),
                       dom.idx=c(rep(NA, 1 + fixed.M + M), 1:dom.M, rep(NA, random.M + random2.M)),
                       random.idx=c(rep(NA, 1 + fixed.M + M + dom.M), 1:random.M, rep(NA, random2.M)),
                       random2.idx=c(rep(NA, 1 + fixed.M + M + dom.M + random.M), 1:random2.M))
        } else if (H == 2) {
          y = c(scale(Y_[1, ], scale=sd(Y_[1, ], na.rm=T)), scale(Y_[2, ], scale=sd(Y_[2, ], na.rm=T)))
          data <- list(y=y,
                       intercept=c(1, rep(NA, fixed.M + M + dom.M + random.M), 2, rep(NA, fixed.M + M + dom.M + random.M)),
                       fixed=c(NA, 1:fixed.M, rep(NA,  M + dom.M + random.M), NA, fixed.M + 1:fixed.M, rep(NA,  M + dom.M + random.M)),
                       idx1=c(rep(NA, 1 +  fixed.M), 1:M, rep(NA, dom.M + random.M), rep(NA, 1 +  fixed.M + M + dom.M + random.M)),
                       idx2=c(rep(NA, 1 +  fixed.M + M + dom.M + random.M), rep(NA, 1 +  fixed.M), 1:M, rep(NA, dom.M + random.M)),
                       dom.idx1=c(rep(NA, 1 + fixed.M + M), 1:dom.M, rep(NA, random.M), rep(NA, 1 + fixed.M + M + dom.M + random.M)),
                       dom.idx2=c(rep(NA, 1 + fixed.M + M + dom.M + random.M), rep(NA, 1 + fixed.M + M),  1:dom.M, rep(NA, random.M)),
                       random.idx=c(rep(NA, 1 + fixed.M + M + dom.M), 1:random.M, rep(NA, 1 + fixed.M + M + dom.M), random.M + 1:random.M))
        }
        formula = NA
        ## http://stackoverflow.com/questions/4951442/formula-with-dynamic-number-of-variables
        pre = "y ~ -1"
        if (fixed.M != 0) {
          fixed = "f(fixed, hyper=list(theta=list(prior=\"loggamma\", param=c(1, 1))))"
          #fixed = "f(fixed, hyper=list(theta=list(prior=\"loggamma\", param=c(1,0.01))))"
          ##fixed = "fixed"
        }

        intercept = "intercept"

        dom.effect = NULL
        random.effect = NULL
        random2.effect = NULL
        if (H == 1) {
          effect = "f(idx1, constr=TRUE, hyper=list(theta=list(prior=\"loggamma\", param=c(1,1))))"
          if (dom.M != 0) {
            dom.effect = "f(dom.idx, hyper=list(theta=list(prior=\"loggamma\", param=c(1,1))))"
          }
        } else {
          effect = "f(idx1, constr=TRUE, hyper=list(theta=list(prior=\"loggamma\", param=c(1,1)))) + f(idx2, constr=TRUE, hyper=list(theta=list(prior=\"loggamma\", param=c(1,1))))"
          if (dom.M != 0) {
            dom.effect = "f(dom.idx1, hyper=list(theta=list(prior=\"loggamma\", param=c(1,1)))) + f(dom.idx2, hyper=list(theta=list(prior=\"loggamma\", param=c(1,1))))"
          }
        }
        if (random.M != 0) {
          if (H == 1) {
            random.effect = "f(random.idx, constr=TRUE, hyper=list(theta=list(prior=\"loggamma\", param=c(1,1))))"
          } else if (H == 2) {
            random.effect = "f(random.idx, constr=TRUE, hyper=list(theta=list(prior=\"loggamma\", param=c(1,1))))"
            ##random.effect = paste("f(random.idx, constr=TRUE, model=\"iid2d\", n=", length(data$intercept),")", sep="")
          }
        }
        if (random2.M != 0) {
          if (H == 1) {
            random2.effect = "f(random2.idx, constr=TRUE, hyper=list(theta=list(prior=\"loggamma\", param=c(1,1))))"
          } else if (H == 2) {
            print("ignore the second random effects")
          }
        }
        all = c(pre, intercept, fixed, effect, dom.effect, random.effect, random2.effect)
        all = all[!is.na(all)]
        formula = as.formula(paste(all, collapse = "+"))
        result = inla(formula,
          data=data, family=family, control.predictor=list(A=x),
          quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975),
          num.threads=num.threads, verbose=F)
        if (length(result$mlik) == 0) {
          print("failed, reattempt..")
          next;
        }
        results$logmliks[case] = result$mlik[1]
        print(result$mlik[1])
        for (i in 1:H) {
          results$betas[case, M_* (i - 1) + 1:M] = result$summary.random[[paste("idx", i, sep="")]]$mean[1:M_]
        }
        results$summary.random[[case]] = result$summary.random
        results$inla.objects[[case]] = list(marginals.random = result$marginals.random, summary.hyperpar=result$summary.hyperpar)
        results$hyper.samples[[case]] = inla.hyperpar.sample(n=1000, result=result, intern=FALSE)

        if (case == S) break;
        case <- case + 1;
      }
      list(results=results, diplotypes=diplotypes)
    },
    get.extra = function(result, H=1, h=1) {
      return (result[[h]]$logmliks)
    },
    is.straineff = function() {
      'Is strain effects model?'
      return (T)
    },
    get.sim.straineff = function(M, result, H=1, h=1) {
      ret = get.straineff(M, result$strain.effects.result)
      ret
    },
    get.samples.mean = function(samples, fieldname, rowidx, colname, mliks) {
      S = length(samples)
      stopifnot(S == length(mliks))
      s = rep(0, S)

      num.effect = length(samples[[i]][[fieldname]][, 'mean'])
      for (i in 1:S) {
        shift = sum(samples[[i]][[fieldname]][, 'mean']) / num.effect
        s[i] = samples[[i]][[fieldname]][rowidx, colname] - shift
      }
      if (use.weight()) {
        return(weighted.mean(s, mliks))
      }
      else {
        return(mean(result$betas[, i]))
      }
    },
    get.ci.from.inla.marginal = function(marginals) {
      effect.num = length(marginals)
      med = rep(0, effect.num)
      mu = rep(0, effect.num)
      hpd.wide.left = rep(0, effect.num)
      hpd.wide.right = rep(0, effect.num)
      hpd.narrow.left = rep(0, effect.num)
      hpd.narrow.right = rep(0, effect.num)
      identity <- function(x) {x}
      for (i in 1:effect.num) {
        med[i] = inla.qmarginal(0.5, marginals[[i]])
        mu[i] = inla.emarginal(identity, marginals[[i]])
        hpd.wide.left[i] = inla.qmarginal(0.025, marginals[[i]])
        hpd.wide.right[i] = inla.qmarginal(0.975, marginals[[i]])
        hpd.narrow.left[i] = inla.qmarginal(0.25, marginals[[i]])
        hpd.narrow.right[i] = inla.qmarginal(0.75, marginals[[i]])
      }
      return(list(med=med, mu=mu,
                  hpd.narrow=cbind(hpd.narrow.left, hpd.narrow.right),
                  hpd.wide=cbind(hpd.wide.left, hpd.wide.right)
                  )
             )
    },
    get.straineff.ci = function(M, result, H=1, h=1) {
      logmliks = as.vector(scale(result$results$logmliks, scale=F))
      mliks = exp(logmliks)
      print (mliks)
      marginal = list()
      for (i in 1:M) {
        marginal[[i]] =
          get.combined.marginal(result$results$inla.objects, paste('idx', h, sep=""),
                                i, mliks)
      }
      return(get.ci.from.inla.marginal(marginal))
    },
    get.deviated.effects = function(M, result, H=1, h=1) {
      get.deviated.effects.mean(M, result)
    },
    get.deviated.effects.mean = function(M, result, H=1, h=1) {
      logmliks = as.vector(scale(result$results$logmliks, scale=F))
      mliks = exp(logmliks)
      rowSums(sapply(1:length(mliks), function(i) (result$results$summary.random[[i]]$dom.idx$mean * mliks[i]))) / sum(mliks)
    },
    get.straineff = function(M, result, H=1, h=1) {
      logmliks = as.vector(scale(result$results$logmliks, scale=F))
      if (use.weight()) {
        mliks = exp(logmliks)
        rowSums(sapply(1:length(mliks), function(i) (result$results$summary.random[[i]]$idx1$mean * mliks[i]))) / sum(mliks)
      }
      else {
        mliks = rep(1, length(logmliks))
        rowSums(sapply(1:length(mliks), function(i) (result$results$summary.random[[i]]$idx1$mean * mliks[i]))) / sum(mliks)
      }
    }
  )
)
