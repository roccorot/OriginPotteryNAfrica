# Load Libraries and setup data ----
library(here)
library(nimble)
library(nimbleCarbon)
library(parallel)
load(here('runscripts','binom','readyrun_binom.RData'))

# Define Parameters ----
constants$distance.a <- d$dist.e.scaled
constants$distance.b <- d$dist.a.scaled
niter = 500000
nburnin = 250000
thin = 25

# Run MCMC ----
this_cluster <- makeCluster(4)
chain_output.m4  <- parLapply(cl=this_cluster,X=1:4,fun=runMCMC_spatialDiffusion.2,dat=dat,inits=inits.dualorigin,code=spatialDiffusion.dualorigin,constants=constants,niter=niter,nburnin=nburnin,thin=thin)
stopCluster(this_cluster)

# Run Diagnostics ----
(rhat.m4 <- coda::gelman.diag(chain_output.m4))
if (any(rhat.m4[[1]][,1]>1.01)) {rhat.m4[[1]][which(rhat.m4[[1]][,1]>1.01),1]}


m4 <- nimbleModel(code = spatialDiffusion.dualorigin, data=dat, constants=constants)
Cm4 <- compileNimble(m4)         # calculateWAIC needs compiled model to exist
samples.m4 <- do.call(rbind, chain_output.m4)    # single matrix of samples
waic.m4 <- calculateWAIC(samples.m4, m4)

# Summarise posteriors ----
m4.res <- list()
m4.res$post <- samples.m4[,grep('beta',colnames(samples.m4))]
(m4.res$waic <- waic.m4)
m4.res$rhat <- rhat.m4
m4.res$ess <- coda::effectiveSize(chain_output.m4)

# Store output ----
save(m4.res,file=here('results_images','binom','results_binom_m4.RData'))
