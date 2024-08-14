# Load Libraries and setup data ----
library(here)
library(nimble)
library(nimbleCarbon)
library(parallel)
load(here('runscripts','binom','readyrun_binom.RData'))

# Define Parameters ----
constants$distance <- d$dist.a.scaled
niter = 500000
nburnin = 250000
thin = 25

# Run MCMC ----
this_cluster <- makeCluster(4)
chain_output.m2  <- parLapply(cl=this_cluster,X=1:4,fun=runMCMC_spatialDiffusion.1,dat=dat,inits=inits.singleorigin,code=spatialDiffusion.singleorigin,constants=constants,niter=niter,nburnin=nburnin,thin=thin)
stopCluster(this_cluster)

# Run Diagnostics ----
(rhat.m2 <- coda::gelman.diag(chain_output.m2))
if (any(rhat.m2[[1]][,1]>1.01)) {rhat.m2[[1]][which(rhat.m2[[1]][,1]>1.01),1]}
if (any(rhat.m2[[1]][1:3,1]>1.01)) {rhat.m2[[1]][which(rhat.m2[[1]][1:3,1]>1.01),1]}


m2 <- nimbleModel(code = spatialDiffusion.singleorigin, data=dat, constants=constants)
Cm2 <- compileNimble(m2)         # calculateWAIC needs compiled model to exist
samples.m2 <- do.call(rbind, chain_output.m2)    # single matrix of samples
waic.m2 <- calculateWAIC(samples.m2, m2)

# Summarise posteriors ----
m2.res <- list()
m2.res$post <- samples.m2[,grep('beta',colnames(samples.m2))]
(m2.res$waic <- waic.m2)
m2.res$rhat <- rhat.m2
m2.res$ess <- coda::effectiveSize(chain_output.m2)

# Store output ----
save(m2.res,file=here('results_images','binom','results_binom_m2.RData'))
