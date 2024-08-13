# Load Libraries and setup data ----
library(here)
library(nimble)
library(nimbleCarbon)
library(parallel)
load(here('runscripts','binom','readyrun_binom.RData'))

# Define Parameters ----
constants$distance <- d$dist.e.scaled
niter = 500000
nburnin = 250000
thin = 25

# Run MCMC ----
this_cluster <- makeCluster(4)
chain_output.m1  <- parLapply(cl=this_cluster,X=1:4,fun=runMCMC_spatialDiffusion.1,dat=dat,inits=inits.singleorigin,code=spatialDiffusion.singleorigin,constants=constants,niter=niter,nburnin=nburnin,thin=thin)
stopCluster(this_cluster)

# Run Diagnostics ----
(rhat.m1 <- coda::gelman.diag(chain_output.m1))
if (any(rhat.m1[[1]][,1]>1.01)) {rhat.m1[[1]][which(rhat.m1[[1]][,1]>1.01),1]}


m1 <- nimbleModel(code = spatialDiffusion.singleorigin, data=dat, constants=constants)
Cm1 <- compileNimble(m1)         # calculateWAIC needs compiled model to exist
samples.m1 <- do.call(rbind, chain_output.m1)    # single matrix of samples
waic.m1 <- calculateWAIC(samples.m1, m1)

# Summarise posteriors ----
m1.res <- list()
m1.res$post <- samples.m1[,grep('beta',colnames(samples.m1))]
(m1.res$waic <- waic.m1)
m1.res$rhat <- rhat.m1
m1.res$ess <- coda::effectiveSize(chain_output.m1)

# Store output ----
save(m1.res,file=here('results_images','binom','results_binom_m1.RData'))
