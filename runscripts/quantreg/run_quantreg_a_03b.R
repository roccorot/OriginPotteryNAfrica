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
niter  <- 300000
nburnin <- 150000
thin <- 15

# Run MCMC ----
# this_cluster <- makeCluster(4)
# chain_output.quant.a  <- parLapply(cl=this_cluster,X=1:4,fun=runMCMCquantreg,dat=dat,inits=inits,code=quantreg,constants=constants,niter=niter,nburnin=nburnin,thin=thin)
# stopCluster(this_cluster)
out.1 <- runMCMCquantreg(seed=1,dat=dat,code=quantreg,constants=constants,inits=inits,niter=niter,nburnin=nburnin,thin=thin)
out.2 <- runMCMCquantreg(seed=2,dat=dat,code=quantreg,constants=constants,inits=inits,niter=niter,nburnin=nburnin,thin=thin)
out.3 <- runMCMCquantreg(seed=3,dat=dat,code=quantreg,constants=constants,inits=inits,niter=niter,nburnin=nburnin,thin=thin)
out.4 <- runMCMCquantreg(seed=4,dat=dat,code=quantreg,constants=constants,inits=inits,niter=niter,nburnin=nburnin,thin=thin)
chaint_output.quant.a <- list(out.1,out.2,out.3,out.4)

# Run Diagnostics ----
(rhat.quant.a <- gelman.diag(chain_output.quant.a))
if (any(rhat.quant.a[[1]][,1]>1.01)) {rhat.quant.a[[1]][which(rhat.quant.a[[1]][,1]>1.01),1]}

# Extract Posterior and Save ----
samples.quant.a  <- do.call(rbind,chain_output.quant.a)
(HPDinterval(mcmc(samples.quant.a[,'beta'])))
res.quant.a <- list()
res.quant.a$post <- samples.quant.a[,!grepl('theta',colnames(samples.quant.a))]
res.quant.a$rhat <- rhat.quant.a
res.quant.a$ess <- effectiveSize(chain_output.quant.a)
save(res.quant.a,file=here('results_images','quantreg','results_quantreg_a.RData'))
