# Load Libraries and setup data ----
library(here)
library(nimble)
library(nimbleCarbon)
library(parallel)
load(here('runscripts','binom','readyrun_binom.RData'))

# Define Parameters ----
constants$distance <- d$dist.o.scaled
niter = 500000
nburnin = 250000
thin = 25

# Run MCMC ----
this_cluster <- makeCluster(4)
chain_output.m3  <- parLapply(cl=this_cluster,X=1:4,fun=runMCMC_spatialDiffusion.1,dat=dat,inits=inits.singleorigin,code=spatialDiffusion.singleorigin,constants=constants,niter=niter,nburnin=nburnin,thin=thin)
stopCluster(this_cluster) 

# Run Diagnostics ----
(rhat.m3 <- coda::gelman.diag(chain_output.m3))
if (any(rhat.m3[[1]][,1]>1.01)) {rhat.m3[[1]][which(rhat.m3[[1]][,1]>1.01),1]}


m3 <- nimbleModel(code = spatialDiffusion.singleorigin, data=dat, constants=constants)
Cm3 <- compileNimble(m3)         # calculateWAIC needs compiled model to exist
samples.m3 <- do.call(rbind, chain_output.m3)    # single matrix of samples
waic.m3 <- calculateWAIC(samples.m3, m3)

# Summarise posteriors ----
m3.res <- list()
m3.res$post <- samples.m3[,grep('beta',colnames(samples.m3))]
(m3.res$waic <- waic.m3)
m3.res$rhat <- rhat.m3
m3.res$ess <- coda::effectiveSize(chain_output.m3)

# Store output ----
save(m3.res,file=here('results_images','binom','results_binom_m3.RData'))
