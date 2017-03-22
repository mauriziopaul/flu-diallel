# 1) Install INLA
# 2) Move contents of Diploffect_experimental/src here
# 3) Install bagpipe backend with `R CMD install --clean --no-docs bagpipe.backend_0.34.tar.gz`
# 4) Updated versions available here: http://valdarlab.unc.edu/software/bagpipe/_build/html/bagpipe.html#r-bagpipe-backend
# 5) To run this script without error, I needed to make the following modification to the file
##   `/usr/local/lib/R/3.3/site-library/INLA/bin/mac/64bit/inla.run`, since the default /usr/lib version of liblzma.5.dylib
##   is the wrong version. Change the line:
##      export DYLD_LIBRARY_PATH="$RHOME_LIB:$DIR:/usr/lib:/opt/local/lib:$DYLD_LIBRARY_PATH"
##   to 
##      export DYLD_LIBRARY_PATH="$RHOME_LIB:$DIR:/usr/local/opt/xz/lib:/usr/lib:/opt/local/lib:$DYLD_LIBRARY_PATH"

source("source.scripts.diploffect.R")
dir.create("plots", recursive=TRUE, showWarnings=FALSE)

####################################### PreCC data
precc.data <- read.csv("Flu-preCC-data.csv", header=TRUE)

precc.D4.inla.diploffect.JAX00072951.model <- run.diploffect.inla.through.genomecache(formula=D4.weight ~ 1 + D0_weight_grams + locus.full, data=precc.data, K=NULL,
                                                                                 genomecache="./preCCflu_happy_cache/", locus="JAX00072951",
                                                                                 num.draws=100, 
                                                                                 seed=1, gamma.rate=1, impute.on="SUBJECT.NAME",
                                                                                 do.augment.of.cache=FALSE)
precc.D4.inla.diploffect.JAX00072951.model.summary <- run.diploffect.inla.summary.stats(precc.D4.inla.diploffect.JAX00072951.model)

pdf("plots/diploffect-preCC.pdf", height=15, width=15)
par(mfrow=c(2,2))
plot.straineff.ci(precc.D4.inla.diploffect.JAX00072951.model.summary, flip=FALSE)
plot.diplotype.ci(precc.D4.inla.diploffect.JAX00072951.model.summary, flip=FALSE)
plot.deviation.ci(precc.D4.inla.diploffect.JAX00072951.model.summary, flip=FALSE)
plot.varexp.ci(precc.D4.inla.diploffect.JAX00072951.model.summary, add.numbers=TRUE)
dev.off()
