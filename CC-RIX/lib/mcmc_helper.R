ifow=function(test, yes, no){ if (test){return (yes)}; no }

mcmc.subset <- function(coda.object, burnin=0, thin=1)
{
	total <- nrow(coda.object[[1]])
	num.wanted <- floor((total-burnin)/thin)
	wanted <- burnin + (1:num.wanted)*thin

	for (i in 1:length(coda.object))
	{
	coda.object[[i]] <- as.mcmc(coda.object[[i]][wanted,])
	}
	as.mcmc.list(coda.object)
}

mcmc.stack <- function(coda.object)
{
	if (inherits(coda.object, "mcmc"))
	{
		return (coda.object)
	}
	if (!inherits(coda.object, "mcmc.list"))
	{
		stop("Non-mcmc object passed to function\n")
	}
	chain <- coda.object[[1]]
	for (i in 2:nchain(coda.object))
	{
		chain <- rbind(chain, coda.object[[i]])
	}
	as.mcmc(chain)
}

plot.hpd <- function(coda.object,
		wanted=varnames(coda.object),
		prob.wide=0.95,
		prob.narrow=0.50,
		xlab="HPD interval",
		names=NULL,
		type="p",
        centered=F,
        pos=0,  ## adjust the label's position
		...)
{
	which.wanted=ifow(is.integer(wanted), wanted, match(wanted, varnames(coda.object)))
	num.wanted=length(which.wanted)

	chain <- mcmc.stack(coda.object)
	mu    <- colMeans(chain[,which.wanted])
	med   <- apply(coda::HPDinterval(chain, prob=0.01)[which.wanted,],
			1, mean)
	hpd.wide    <- coda::HPDinterval(chain, prob=prob.wide)[which.wanted,]
	hpd.narrow  <- coda::HPDinterval(chain, prob=prob.narrow)[which.wanted,]

	if (is.null(names)) names <- varnames(chain)[which.wanted]
	else names <- rep(names, length.out=length(wanted))
	ypos <- plot.ci(med, hpd.narrow, hpd.wide, names=names, xlab=xlab, col.midvals="white", pch.midvals="|", type=type, pos=pos, ...)
	if ("p"==type)
	{
		points(mu, ypos, pch="|")
	}
	invisible(ypos)
}
