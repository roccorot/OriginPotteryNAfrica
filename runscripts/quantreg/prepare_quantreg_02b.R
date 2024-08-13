# Load Libraries and Data ----
library(here)
library(nimble)
library(nimbleCarbon);data(intcal20)
library(coda)
load(here('data','input_quantreg.RData'))

# Setup lists and model code for nimble ----
# Data
dat  <- list()
dat$cra  <- d$Yrs.BP
dat$cra.error <- d$STD
dat$chronological.constraint <- 1

# Constants
constants <- list()
constants$calBP  <- intcal20$CalBP
constants$C14BP  <- intcal20$C14Age
constants$C14err <- intcal20$C14Age.sigma
constants$n <- nrow(d)
constants$tau <- 0.95

# Core Code 
quantreg  <- nimbleCode({
	for (i in 1:n)
	{
	mu[i] <- alpha + beta*distance[i]
	theta[i] ~ dAsymLaplace(mu=mu[i],sigma=sigma,tau=tau)
	c14age[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
	sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
	sigmaDate[i] <- (cra.error[i]^2+sigmaCurve[i]^2)^(1/2);
	cra[i] ~ dnorm(mean=c14age[i],sd=sigmaDate[i]);
	}
	alpha ~ T(dnorm(11000,sd=1000),500,50000)
	beta ~ dnorm(0,sd=1)
	sigma ~ dexp(0.05) #should be smaller?
	constraint.text.to.replace
})

quantreg <- gsub('constraint.text.to.replace',constraint.text,deparse(quantreg)) |> parse(text=_)

# Inits 
inits  <- list()
inits$alpha <- 11000
inits$beta <- 0
inits$sigma <- 100
inits$theta <- d$theta

# Create run function for executing MCMC in parallel ----
runMCMCquantreg  <- function(seed,dat,code,constants,inits,niter,nburnin,thin)
{
	library(nimbleCarbon)
	fit.model  <- nimbleModel(code=code,constants=constants,data=dat,inits=inits)
	cfit.model  <- compileNimble(fit.model)
	MCMC <- buildMCMC(cfit.model,monitors=c('alpha','beta','sigma','theta'))
	cMCMC <- compileNimble(MCMC)
	results <- runMCMC(cMCMC,niter=niter,thin=thin,nburnin=nburnin,setSeed=seed,samplesAsCodaMCMC=TRUE)
}

# Save relevant objects into an R image file ----
save(runMCMCquantreg,quantreg,constants,dat,inits,file=here('runscripts','quantreg','readyrun_quantreg.RData'))
