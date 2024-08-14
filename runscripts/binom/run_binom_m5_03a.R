# Load Libraries and setup data ----
library(here)
library(nimble)
library(nimbleCarbon)
library(parallel)
load(here('runscripts','binom','readyrun_binom.RData'))

# Define Parameters ----
constants$distance.a <- d$dist.e.scaled
constants$distance.b <- d$dist.o.scaled
niter = 500000
nburnin = 250000
thin = 25

# Run MCMC ----
this_cluster <- makeCluster(4)
chain_output.m5  <- parLapply(cl=this_cluster,X=1:4,fun=runMCMC_spatialDiffusion.2,dat=dat,inits=inits.dualorigin,code=spatialDiffusion.dualorigin,constants=constants,niter=niter,nburnin=nburnin,thin=thin)
stopCluster(this_cluster)

# Run Diagnostics ----
(rhat.m5 <- coda::gelman.diag(chain_output.m5))
if (any(rhat.m5[[1]][,1]>1.01)) {rhat.m5[[1]][which(rhat.m5[[1]][,1]>1.01),1]}
if (any(rhat.m5[[1]][1:6,1]>1.01)) {rhat.m5[[1]][which(rhat.m5[[1]][1:6,1]>1.01),1]}


m5 <- nimbleModel(code = spatialDiffusion.dualorigin, data=dat, constants=constants)
Cm5 <- compileNimble(m5)         # calculateWAIC needs compiled model to exist
samples.m5 <- do.call(rbind, chain_output.m5)    # single matrix of samples
waic.m5 <- calculateWAIC(samples.m5, m5)

# Summarise posteriors ----
m5.res <- list()
m5.res$post <- samples.m5[,grep('beta',colnames(samples.m5))]
(m5.res$waic <- waic.m5)
m5.res$rhat <- rhat.m5
m5.res$ess <- coda::effectiveSize(chain_output.m5)

# Store output ----
save(m5.res,file=here('results_images','binom','results_binom_m5.RData'))
