# Load Libraries and setup data ----
library(nimble)
library(nimbleCarbon)
library(here)
library(coda)
library(parallel)

load(here('runscripts','quantreg','readyrun_quantreg.RData'))
load(here('data','input_quantreg.RData'))

# Define Parameters ----
constants$distance <- d$dist.e
niter  <- 100000
nburnin <- 50000
thin <- 5

# Run MCMC ----
this_cluster <- makeCluster(4)
chain_output.quant.e  <- parLapply(cl=this_cluster,X=1:4,fun=runMCMCquantreg,dat=dat,inits=inits,code=quantreg,constants=constants,niter=niter,nburnin=nburnin,thin=thin)
stopCluster(this_cluster)

# Run Diagnostics ----
(rhat.quant.e <- gelman.diag(chain_output.quant.e))
if (any(rhat.quant.e[[1]][,1]>1.01)) {rhat.quant.e[[1]][which(rhat.quant.e[[1]][,1]>1.01),1]}

# Extract Posterior and Save ----
samples.quant.e  <- do.call(rbind,chain_output.quant.e)
(HPDinterval(mcmc(samples.quant.e[,'beta'])))
res.quant.e <- list()
res.quant.e$post <- samples.quant.e[,!grepl('theta',colnames(samples.quant.e))]
res.quant.e$rhat <- rhat.quant.e
res.quant.e$ess <- effectiveSize(chain_output.quant.e)
save(res.quant.e,file=here('results_images','quantreg','results_quantreg_e.RData'))
