library(here)
library(coda)
library(ggplot2)
library(patchwork)
library(purrr)
library(viridis)

logistic  <- function(x){p <- 1/(1+exp(-x))}


# Load Observed Data and Results ----
load(here('data','input_binom.RData'))
load(here('results_images','binom','results_binom_m1.RData'))
load(here('results_images','binom','results_binom_m2.RData'))
load(here('results_images','binom','results_binom_m3.RData'))
load(here('results_images','binom','results_binom_m4.RData'))
load(here('results_images','binom','results_binom_m5.RData'))
load(here('results_images','binom','results_binom_m6.RData'))
load(here('results_images','binom','results_binom_m7.RData'))

mean.dists <- mean(c(d$dist.e,d$dist.a,d$dist.o))
sd.dists <- sd(c(d$dist.e,d$dist.a,d$dist.o))
mean.theta <- mean(d$theta)
sd.theta <- sd(d$theta)


# Marginal Model Prediction Plots ----
xrange <- seq(0,3000,length.out=50)
yrange <- seq(7000,11000,length.out=50)
pred.mat <- expand.grid(x.unscaled=xrange,y.unscaled=yrange)
pred.mat$x <- ((pred.mat$x.unscaled*1000) - mean.dists)/sd.dists
pred.mat$y <- (pred.mat$y.unscaled-mean.theta)/sd.theta

set.seed(123)
nsim <- 1000
index <- sample(1:nrow(m1.res$post),size=1000)
posts <- list(m1.res$post[index,],m2.res$post[index,],m3.res$post[index,],m4.res$post[index,],m5.res$post[index,],m6.res$post[index,],m7.res$post[index,])
combinations <- cbind(c(T,F,F,T,T,F,T),c(F,T,F,T,F,T,T),c(F,F,T,F,T,T,T))
colnames(combinations) <- c('e','a','o')
rownames(combinations) <- paste0('m',1:7)
model.names <- paste0('m',1:7)
origins <- c('Bir Kiseiba','Adrar Bous','Onjoungo Raven de la Mouch')
model.slots <- vector('list',length=7)
for(i in 1:7){model.slots[[i]] <- vector('list',length=3)}

for (i in 1:7)
{
	for (j in 1:3)
	{
		if (combinations[i,j]==TRUE)
		{
			tmp <- pred.mat
			if (i %in% c(1:3))
			{
				for (k in 1:nrow(tmp))
				{
				tmp$p[k] <- apply(posts[[i]],1,function(x,y,z){return(logistic(x[1]-x[2]*y-x[3]*z))},y=tmp$y[k],z=tmp$x[k]) |> median()
				}

			}
			if (i %in% c(4:6))
			{
				if(j==which(combinations[i,])[1])
				{
					for (k in 1:nrow(tmp))
					{
						tmp$p[k] <- apply(posts[[i]],1,function(x,y,z){return(logistic(x[1]-x[3]*y-x[5]*z))},y=tmp$y[k],z=tmp$x[k]) |> median()
					}
				}

				if(j==which(combinations[i,])[2])
				{
					for (k in 1:nrow(tmp))
					{
						tmp$p[k] <- apply(posts[[i]],1,function(x,y,z){return(logistic(x[2]-x[4]*y-x[6]*z))},y=tmp$y[k],z=tmp$x[k]) |> median()
					}
				}
			}

			if (i==7)
			{
				if (j==1)
				{
					for (k in 1:nrow(tmp))
					{
						tmp$p[k] <- apply(posts[[i]],1,function(x,y,z){return(logistic(x[1]-x[4]*y-x[7]*z))},y=tmp$y[k],z=tmp$x[k]) |> median()
					}

				}

				if (j==2)
				{
					for (k in 1:nrow(tmp))
					{
						tmp$p[k] <- apply(posts[[i]],1,function(x,y,z){return(logistic(x[2]-x[5]*y-x[8]*z))},y=tmp$y[k],z=tmp$x[k]) |> median()
					}

				}

				if (j==3)
				{
					for (k in 1:nrow(tmp))
					{
						tmp$p[k] <- apply(posts[[i]],1,function(x,y,z){return(logistic(x[3]-x[6]*y-x[9]*z))},y=tmp$y[k],z=tmp$x[k]) |> median()
					}
				}
			}
			model.slots[[i]][[j]]  <- ggplot() + 
				geom_tile(data=tmp,aes(x=x.unscaled,y=y.unscaled,fill=p)) +
				scale_fill_viridis(limits=c(0,1)) +
				labs(title = paste0(model.names[i],'/',origins[j]),x='Distance from Origin (km)',y='Cal BP') +
				scale_y_reverse() +
				theme(legend.position='bottom') +
				theme_minimal()
		}

		if (combinations[i,j]==FALSE)
		{
			model.slots[[i]][[j]] <- ggplot() + geom_blank() + theme_minimal()

		}
	}
}

model.slots.flattened <- flatten(model.slots)
combined <- reduce(model.slots.flattened,`+`)
pdf(width=10,height=16,file=here('figures','si','marginal_predictions.pdf'))
combined + plot_layout(guides='collect',ncol=3)
dev.off()

# Posterior Summaries ----
post.summary <- data.frame(models = rep(paste0('m',1:7),c(3,3,3,6,6,6,9)))
post.summary$param <- c('beta.e.0','beta.e.time','beta.e.distance',
			'beta.a.0','beta.a.time','beta.a.distance',
			'beta.o.0','beta.o.time','beta.o.time',
			'beta.e.0','beta.e.time','beta.e.distance','beta.a.0','beta.a.time','beta.a.distance',
			'beta.e.0','beta.e.time','beta.e.distance','beta.o.0','beta.o.time','beta.o.distance',
			'beta.a.0','beta.a.time','beta.a.distance','beta.o.0','beta.o.time','beta.o.distance',
			'beta.e.0','beta.e.time','beta.e.distance','beta.a.0','beta.a.time','beta.a.distance','beta.o.0','beta.o.time','beta.o.distance')

posts <- cbind(m1.res$post,
	       m2.res$post,
	       m3.res$post,
	       m4.res$post[,c(1,3,5,2,4,6)],
	       m5.res$post[,c(1,3,5,2,4,6)],
	       m6.res$post[,c(1,3,5,2,4,6)],
	       m7.res$post[,c(1,4,7,2,5,8,3,6,9)])
rhats <- c(m1.res$rhat[[1]][1:3,1],
	   m2.res$rhat[[1]][1:3,1],
	   m3.res$rhat[[1]][1:3,1],
	   m4.res$rhat[[1]][c(1,3,5,2,4,6),1],
	   m5.res$rhat[[1]][c(1,3,5,2,4,6),1],
	   m6.res$rhat[[1]][c(1,3,5,2,4,6),1],
	   m6.res$rhat[[1]][c(1,4,7,2,5,8,3,6,9),1])

ess <- c(m1.res$ess[1:3],
	  m2.res$ess[1:3],
	  m3.res$ess[1:3],
	  m4.res$ess[c(1,3,5,2,4,6)],
	  m5.res$ess[c(1,3,5,2,4,6)],
	  m6.res$ess[c(1,3,5,2,4,6)],
	  m7.res$ess[c(1,4,7,2,5,8,3,6,9)])

post.summary$ess <- post.summary$rhat  <- post.summary$hpdi  <- post.summary$mean <- NA
param.numbers <- c(3,3,3,6,6,6,9)
model.numbers <- 1:7
current.index <- 1

for (i in model.numbers)
{
	fill.index <- current.index:c(current.index+param.numbers[i]-1)
	current.index <- max(fill.index) + 1
	post.summary$mean[fill.index] <- round(apply(posts[,fill.index],2,mean),5)
	post.summary$hpdi[fill.index] <- apply(posts[,fill.index],2,function(x){paste0(round(as.numeric(HPDinterval(mcmc(x))),5),collapse=' ~ ')})
	post.summary$rhat[fill.index] <- round(rhats[fill.index],4)
	post.summary$ess[fill.index] <- round(ess[fill.index])
}

write.csv(post.summary,file=here('tables','si','binomial_posterior.csv'))
