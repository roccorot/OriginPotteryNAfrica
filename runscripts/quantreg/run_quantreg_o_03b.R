# Load Libraries and setup data ----
library(nimble)
library(nimbleCarbon)
library(here)
library(coda)
library(parallel)

load(here('runscripts','quantreg','readyrun_quantreg.RData'))
load(here('data','input_quantreg.RData'))

# Define Parameters
constants$distance <- d$dist.o
niter  <- 500000
nburnin <- 250000
thin <- 25

# Run MCMC ----
this_cluster <- makeCluster(4)
chain_output.quant.o  <- parLapply(cl=this_cluster,X=1:4,fun=runMCMCquantreg,dat=dat,inits=inits,code=quantreg,constants=constants,niter=niter,nburnin=nburnin,thin=thin)
stopCluster(this_cluster)

# Run Diagnostics ----
(rhat.quant.o <- gelman.diag(chain_output.quant.o))
if (any(rhat.quant.o[[1]][,1]>1.01)) {rhat.quant.o[[1]][which(rhat.quant.o[[1]][,1]>1.01),1]}

# Extract Posterior and Save ----
samples.quant.o  <- do.call(rbind,chain_output.quant.o)
(HPDinterval(mcmc(samples.quant.o[,'beta'])))
res.quant.o <- list()
res.quant.o$post <- samples.quant.o[,!grepl('theta',colnames(samples.quant.o))]
res.quant.o$rhat <- rhat.quant.o
res.quant.o$ess <- effectiveSize(chain_output.quant.o)
save(res.quant.o,file=here('results_images','quantreg','results_quantreg_o.RData'))
