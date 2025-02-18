# Load Libraries and setup data ----
library(nimble)
library(nimbleCarbon)
library(here)
library(coda)
library(parallel)

load(here('runscripts','quantreg','readyrun_quantreg.RData'))
load(here('data','input_quantreg.RData'))

# Define Parameters ----
constants$distance <- d$dist.a
niter  <- 100000
nburnin <- 50000
thin <- 5

# Run MCMC ----
this_cluster <- makeCluster(4)
chain_output.quant.a  <- parLapply(cl=this_cluster,X=1:4,fun=runMCMCquantreg,dat=dat,inits=inits,code=quantreg,constants=constants,niter=niter,nburnin=nburnin,thin=thin)
stopCluster(this_cluster)

# Run Diagnostics ----
(rhat.quant.a <- gelman.diag(chain_output.quant.a))
if (any(rhat.quant.a[[1]][,1]>1.01)) {rhat.quant.a[[1]][which(rhat.quant.a[[1]][,1]>1.01),1]}
# Rhat>1.01 for a theta112, can be neglected as 14C dates may not have a gaussian posterior

# Extract Posterior and Save ----
samples.quant.a  <- do.call(rbind,chain_output.quant.a)
(HPDinterval(mcmc(samples.quant.a[,'beta'])))
res.quant.a <- list()
res.quant.a$post <- samples.quant.a[,!grepl('theta',colnames(samples.quant.a))]
res.quant.a$rhat <- rhat.quant.a
res.quant.a$ess <- effectiveSize(chain_output.quant.a)
save(res.quant.a,file=here('results_images','quantreg','results_quantreg_a.RData'))
