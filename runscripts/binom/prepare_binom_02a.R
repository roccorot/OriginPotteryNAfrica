# Load library and data ----
library(here)
library(dplyr)
library(rcarbon)
library(nimbleCarbon);data(intcal20)
load(here('data','input_binom.RData'))

# Define constants ----
# constants 
constants <- list()
constants$n.dates <- nrow(d)
constants$calBP  <- intcal20$CalBP
constants$C14BP  <- intcal20$C14Age
constants$C14err <- intcal20$C14Age.sigma
constants$m.theta  <- mean(d$theta)
constants$sd.theta <- sd(d$theta)

# Define Data ----
dat <- list()
dat$cra <- d$Yrs.BP
dat$cra.error <- d$STD
dat$pottery <- d$pottery_bin 
dat$chronological.constraint <- 1

# Define Spatio-Temporal Diffusion Models ----
spatialDiffusion.singleorigin <- nimbleCode({
	for (i in 1:n.dates)
	{
		pottery[i] ~ dbern(p[i]);
		logit(p[i])  <- beta0 - beta1 * scaled.theta[i] - beta2 * distance[i];
		scaled.theta[i] <- (theta[i] - m.theta)/sd.theta
		theta[i] ~ dunif(1000,30000);
		eta[i]  <- interpLin(z=theta[i],x=calBP[],y=C14BP[]);
		sigmaCurve[i] <- interpLin(z=theta[i],x=calBP[],y=C14err[]);
		sigma[i]  <- (cra.error[i]^2+sigmaCurve[i]^2)^(1/2);
		cra[i] ~ dnorm(mean=eta[i],sd=sigma[i]);
	}
	beta0 ~ dnorm(0,sd=0.5); #there are probably better priors
	beta1 ~ T(dnorm(mean=0,sd=0.5),0,Inf)
	beta2 ~ T(dnorm(mean=0,sd=0.5),0,Inf)
	constraint.text.to.replace
})
spatialDiffusion.singleorigin <- gsub('constraint.text.to.replace',constraint.text,deparse(spatialDiffusion.singleorigin)) |> parse(text=_)

spatialDiffusion.dualorigin <- nimbleCode({
	for (i in 1:n.dates)
	{
		pottery[i] ~ dbern(max(p[i,1:2]));
		logit(p[i,1])  <- beta0.a - beta1.a * scaled.theta[i] - beta2.a * distance.a[i];
		logit(p[i,2])  <- beta0.b - beta1.b * scaled.theta[i] - beta2.b * distance.b[i];
		scaled.theta[i] <- (theta[i] - m.theta)/sd.theta
		theta[i] ~ dunif(1000,30000);
		eta[i]  <- interpLin(z=theta[i],x=calBP[],y=C14BP[]);
		sigmaCurve[i] <- interpLin(z=theta[i],x=calBP[],y=C14err[]);
		sigma[i]  <- (cra.error[i]^2+sigmaCurve[i]^2)^(1/2);
		cra[i] ~ dnorm(mean=eta[i],sd=sigma[i]);
	}

	beta0.a ~ dnorm(0,sd=0.5);
	beta0.b ~ dnorm(0,sd=0.5); #there are probably better priors
	beta1.a ~ T(dnorm(mean=0,sd=0.5),0,Inf)
	beta1.b ~ T(dnorm(mean=0,sd=0.5),0,Inf)
	beta2.a ~ T(dnorm(mean=0,sd=0.5),0,Inf)
	beta2.b ~ T(dnorm(mean=0,sd=0.5),0,Inf)
	constraint.text.to.replace
})
spatialDiffusion.dualorigin <- gsub('constraint.text.to.replace',constraint.text,deparse(spatialDiffusion.dualorigin)) |> parse(text=_)

spatialDiffusion.tripleorigin <- nimbleCode({
	for (i in 1:n.dates)
	{
		pottery[i] ~ dbern(max(p[i,1:3]));
		logit(p[i,1])  <- beta0.a - beta1.a * scaled.theta[i] - beta2.a * distance.a[i];
		logit(p[i,2])  <- beta0.b - beta1.b * scaled.theta[i] - beta2.b * distance.b[i];
		logit(p[i,3])  <- beta0.c - beta1.c * scaled.theta[i] - beta2.c * distance.c[i];
		scaled.theta[i] <- (theta[i] - m.theta)/sd.theta
		theta[i] ~ dunif(1000,30000);
		eta[i]  <- interpLin(z=theta[i],x=calBP[],y=C14BP[]);
		sigmaCurve[i] <- interpLin(z=theta[i],x=calBP[],y=C14err[]);
		sigma[i]  <- (cra.error[i]^2+sigmaCurve[i]^2)^(1/2);
		cra[i] ~ dnorm(mean=eta[i],sd=sigma[i]);
	}

	beta0.a ~ dnorm(0,sd=0.5);
	beta0.b ~ dnorm(0,sd=0.5); #there are probably better priors
	beta0.c ~ dnorm(0,sd=0.5); #there are probably better priors
	beta1.a ~ T(dnorm(mean=0,sd=0.5),0,Inf) 
	beta1.b ~ T(dnorm(mean=0,sd=0.5),0,Inf) 
	beta1.c ~ T(dnorm(mean=0,sd=0.5),0,Inf) 
	beta2.a ~ T(dnorm(mean=0,sd=0.5),0,Inf) 
	beta2.b ~ T(dnorm(mean=0,sd=0.5),0,Inf) 
	beta2.c ~ T(dnorm(mean=0,sd=0.5),0,Inf) 
	constraint.text.to.replace
})
spatialDiffusion.tripleorigin <- gsub('constraint.text.to.replace',constraint.text,deparse(spatialDiffusion.tripleorigin)) |> parse(text=_)



# Define init values ----
inits.singleorigin  <- inits.dualorigin <- inits.tripleorigin <- list()
inits.singleorigin$beta0 <- 0
inits.singleorigin$beta1 <- 0
inits.singleorigin$beta2 <- 0
inits.singleorigin$theta <- d$theta

inits.dualorigin$beta0.a <- inits.dualorigin$beta0.b <- 0
inits.dualorigin$beta1.a <- inits.dualorigin$beta1.b <- 0
inits.dualorigin$beta2.a <- inits.dualorigin$beta2.b <- 0
inits.dualorigin$theta <- d$theta

inits.tripleorigin$beta0.a <- inits.tripleorigin$beta0.b <- inits.tripleorigin$beta0.c <- 0
inits.tripleorigin$beta1.a <- inits.tripleorigin$beta1.b <- inits.tripleorigin$beta1.c <- 0
inits.tripleorigin$beta2.a <- inits.tripleorigin$beta2.b <- inits.tripleorigin$beta2.c <- 0
inits.tripleorigin$theta <- d$theta

# Generate Run functions  ----
runMCMC_spatialDiffusion.1  <- function(seed,dat,code,constants,inits,useWAIC=TRUE,niter,nburnin,thin)
{
	library(nimbleCarbon)
	fit.model  <- nimbleModel(code=code,constants=constants,data=dat,inits=inits)
	cfit.model  <- compileNimble(fit.model)
	if (useWAIC)
	{
		monitors <- fit.model$getParents(fit.model$getNodeNames(dataOnly = TRUE), stochOnly = TRUE)
	}

	MCMC <- buildMCMC(cfit.model,monitors=monitors)
	cMCMC <- compileNimble(MCMC)
	results <- runMCMC(cMCMC,niter=niter,thin=thin,nburnin=nburnin,setSeed=seed,samplesAsCodaMCMC=TRUE)
}


runMCMC_spatialDiffusion.2  <- function(seed,dat,code,constants,inits,useWAIC=TRUE,niter,nburnin,thin)
{
	library(nimbleCarbon)
	fit.model  <- nimbleModel(code=code,constants=constants,data=dat,inits=inits)
	cfit.model  <- compileNimble(fit.model)
	if (useWAIC)
	{
		monitors <- fit.model$getParents(fit.model$getNodeNames(dataOnly = TRUE), stochOnly = TRUE)
	}
	conf <- configureMCMC(cfit.model,monitors=monitors)
	conf$removeSampler(c('beta0.a','beta0.b','beta1.a','beta1.b','beta2.a','beta2.b'))
	conf$addSampler(c('beta0.a','beta0.b','beta1.a','beta1.b','beta2.a','beta2.b'),type='AF_slice',control=list(sliceAdaptFactorInterval=500,sliceAdaptWidthTolerance=0.1,sliceAdaptFactorMaxIter=30000))

	MCMC <- buildMCMC(conf)
	cMCMC <- compileNimble(MCMC)
	results <- runMCMC(cMCMC,niter=niter,thin=thin,nburnin=nburnin,setSeed=seed,samplesAsCodaMCMC=TRUE)
}


runMCMC_spatialDiffusion.3  <- function(seed,dat,code,constants,inits,useWAIC=TRUE,niter,nburnin,thin)
{
	library(nimbleCarbon)
	fit.model  <- nimbleModel(code=code,constants=constants,data=dat,inits=inits)
	cfit.model  <- compileNimble(fit.model)
	if (useWAIC)
	{
		monitors <- fit.model$getParents(fit.model$getNodeNames(dataOnly = TRUE), stochOnly = TRUE)
	}
	conf <- configureMCMC(cfit.model,monitors=monitors)
	conf$removeSampler(c('beta0.a','beta0.b','beta0.c','beta1.a','beta1.b','beta1.c','beta2.a','beta2.b','beta2.c'))
	conf$addSampler(c('beta0.a','beta0.b','beta0.c','beta1.a','beta1.b','beta1.c','beta2.a','beta2.b','beta2.c'),type='AF_slice',control=list(sliceAdaptFactorInterval=500,sliceAdaptWidthTolerance=0.1,sliceAdaptFactorMaxIter=30000))
	MCMC <- buildMCMC(conf)
	cMCMC <- compileNimble(MCMC)
	results <- runMCMC(cMCMC,niter=niter,thin=thin,nburnin=nburnin,setSeed=seed,samplesAsCodaMCMC=TRUE)
}

# Store everything in an R image file ----
save(constants,dat,d,inits.singleorigin,inits.dualorigin,inits.tripleorigin,spatialDiffusion.singleorigin,spatialDiffusion.dualorigin,spatialDiffusion.tripleorigin,runMCMC_spatialDiffusion.1,runMCMC_spatialDiffusion.2,runMCMC_spatialDiffusion.3,file=here('runscripts','binom','readyrun_binom.RData'))
