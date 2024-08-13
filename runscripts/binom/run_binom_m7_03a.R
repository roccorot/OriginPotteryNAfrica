# Load Libraries and setup data ----
library(here)
library(nimble)
library(nimbleCarbon)
library(parallel)
load(here('runscripts','binom','readyrun_binom.RData'))

# Define Parameters ----
constants$distance.a <- d$dist.e.scaled
constants$distance.b <- d$dist.a.scaled
constants$distance.c <- d$dist.o.scaled
niter = 500000
nburnin = 250000
thin = 25

# Run MCMC ----
this_cluster <- makeCluster(4)
chain_output.m7  <- parLapply(cl=this_cluster,X=1:4,fun=runMCMC_spatialDiffusion.3,dat=dat,inits=inits.tripleorigin,code=spatialDiffusion.tripleorigin,constants=constants,niter=niter,nburnin=nburnin,thin=thin)
stopCluster(this_cluster)

# Run Diagnostics ----
(rhat.m7 <- coda::gelman.diag(chain_output.m7))
if (any(rhat.m7[[1]][,1]>1.01)) {rhat.m7[[1]][which(rhat.m7[[1]][,1]>1.01),1]}


m7 <- nimbleModel(code = spatialDiffusion.tripleorigin, data=dat, constants=constants)
Cm7 <- compileNimble(m7)         # calculateWAIC needs compiled model to exist
samples.m7 <- do.call(rbind, chain_output.m7)    # single matrix of samples
waic.m7 <- calculateWAIC(samples.m7, m7)

# Summarise posteriors ----
m7.res <- list()
m7.res$post <- samples.m7[,grep('beta',colnames(samples.m7))]
(m7.res$waic <- waic.m7)
m7.res$rhat <- rhat.m7
m7.res$ess <- coda::effectiveSize(chain_output.m7)

# Store output ----
save(m7.res,file=here('results_images','binom','results_binom_m7.RData'))
