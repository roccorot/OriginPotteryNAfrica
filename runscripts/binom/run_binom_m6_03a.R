# Load Libraries and setup data ----
library(here)
library(nimble)
library(nimbleCarbon)
library(parallel)
load(here('runscripts','binom','readyrun_binom.RData'))

# Define Parameters ----
constants$distance.a <- d$dist.a.scaled
constants$distance.b <- d$dist.o.scaled
niter = 500000
nburnin = 250000
thin = 25

# Run MCMC ----
this_cluster <- makeCluster(4)
chain_output.m6  <- parLapply(cl=this_cluster,X=1:4,fun=runMCMC_spatialDiffusion.2,dat=dat,inits=inits.dualorigin,code=spatialDiffusion.dualorigin,constants=constants,niter=niter,nburnin=nburnin,thin=thin)
stopCluster(this_cluster)

# Run Diagnostics ----
(rhat.m6 <- coda::gelman.diag(chain_output.m6))
if (any(rhat.m6[[1]][,1]>1.01)) {rhat.m6[[1]][which(rhat.m6[[1]][,1]>1.01),1]}


m6 <- nimbleModel(code = spatialDiffusion.dualorigin, data=dat, constants=constants)
Cm6 <- compileNimble(m6)         # calculateWAIC needs compiled model to exist
samples.m6 <- do.call(rbind, chain_output.m6)    # single matrix of samples
waic.m6 <- calculateWAIC(samples.m6, m6)

# Summarise posteriors ----
m6.res <- list()
m6.res$post <- samples.m6[,grep('beta',colnames(samples.m6))]
(m6.res$waic <- waic.m6)
m6.res$rhat <- rhat.m6
m6.res$ess <- coda::effectiveSize(chain_output.m6)

# Store output ----
save(m6.res,file=here('results_images','binom','results_binom_m6.RData'))
