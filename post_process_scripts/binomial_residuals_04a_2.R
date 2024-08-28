# Load Library, Data, and Results ----
library(here)
library(rcarbon)
library(sf)
library(spdep)
library(sfdep)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(rnaturalearth)

load(here('data','input_binom.RData'))
load(here('results_images','binom','results_binom_m5.RData'))
load(here('results_images','binom','results_binom_m7.RData'))
background <- ne_countries(scale = "medium", returnclass = "sf") |> st_geometry()

# Setup -----
nsim <- 1000
set.seed(123)
index <- sample(1:nrow(m5.res$post),size=nsim)
post.m7 <- m7.res$post[index,]
post.m5 <- m5.res$post[index,]

# Make predictions and compute residuals ----
predictions.m5 <- sapply(1:nsim,function(x,post,t,d1,d2){return(pmax(c(1/(1+exp(-(post[x,1] - post[x,3]*t - post[x,5]*d1)))),c(1/(1+exp(-(post[x,2] - post[x,4]*t - post[x,6]*d2))))))},post=post.m5,t=d$theta.scaled,d1=d$dist.e.scaled,d2=d$dist.o.scaled)

predictions.m7 <- sapply(1:nsim,function(x,post,t,d1,d2,d3){return(pmax(c(1/(1+exp(-(post[x,1] - post[x,4]*t - post[x,7]*d1)))),c(1/(1+exp(-(post[x,2] - post[x,5]*t - post[x,8]*d2)))),c(1/(1+exp(-(post[x,3] - post[x,6]*t - post[x,9]*d3))))))},post=post.m7,t=d$theta.scaled,d1=d$dist.e.scaled,d2=d$dist.a.scaled,d3=d$dist.o.scaled)


d$residuals.m7  <- d$pottery_bin - apply(predictions.m5,1,median)
d$residuals.m5  <- d$pottery_bin - apply(predictions.m7,1,median)
d.calibrated  <- calibrate(d$Yrs.BP,d$STD)

# Plot Residuals  -----
sites <- vector('list',length=3)
sites[[1]] <- st_as_sf(d[which.CalDates(d.calibrated,BP < (9000+500) & BP > (9000 - 500),p=0.5),],coords=c('Longitude','Latitude'),crs=4326)
sites[[2]] <- st_as_sf(d[which.CalDates(d.calibrated,BP < (8000+500) & BP > (8000 - 500),p=0.5),],coords=c('Longitude','Latitude'),crs=4326)
sites[[3]] <- st_as_sf(d[which.CalDates(d.calibrated,BP < (7000+500) & BP > (7000 - 500),p=0.5),],coords=c('Longitude','Latitude'),crs=4326)
time.slices <- c(9000,8000,7000)
ggm5 <- ggm7 <- vector('list',length=3)

for(i in 1:3)
{
	nb <- st_dist_band(sites[[i]],upper=1000)
	wt  <- st_weights(nb,allow_zero=TRUE)
	Gi.m5.8k = local_g_perm(sites[[i]]$residuals.m5, nb, wt, nsim = 999)
	Gi.m7.8k = local_g_perm(sites[[i]]$residuals.m7, nb, wt, nsim = 999)
	sites[[i]]$cluster.m5 <- Gi.m5.8k$cluster
	sites[[i]]$pval.m5 <- ifelse(Gi.m5.8k$p_folded_sim<=0.05,'<.05','>.05')
	sites[[i]]$cluster.m7 <- Gi.m7.8k$cluster
	sites[[i]]$pval.m7 <- ifelse(Gi.m7.8k$p_folded_sim<=0.05,'<.05','>.05')
	ggm5[[i]] <- ggplot() +
		geom_sf(data=background) +
		geom_sf(data = sites[[i]],aes(color=cluster.m5,shape=pval.m5),size=1.5) +
		scale_fill_manual(values=c('high'='red','low'='blue'))+
		scale_shape_manual(values=c('<.05'=19,'>.05'=21),drop=F,limits=c('<.05','>.05'))+
		xlim(-17.5,40.5) +
		ylim(7.5,34.5) +
		labs(title=paste0('m5/',time.slices[[i]],'BP'),color='Gi-cluster',shape='P-value',x='Longitude',y='Latitude')+
		theme_minimal() +
		theme(plot.title = element_text(size=9), axis.title = element_text(size=8),axis.text = element_text(size=8))

	ggm7[[i]] <- ggplot() +
		geom_sf(data=background) +
		geom_sf(data = sites[[i]],aes(color=cluster.m7,shape=pval.m7),size=1.5) +
		scale_fill_manual(values=c('high'='red','low'='blue'))+
		scale_shape_manual(values=c('<.05'=19,'>.05'=21),drop=F,limits=c('<.05','>.05'))+
		xlim(-17.5,40.5) +
		ylim(7.5,34.5) +
		labs(title=paste0('m7/',time.slices[[i]],'BP'),color='Gi-cluster',shape='P-value',x='Longitude',y='Latitude')+
		theme_minimal() +
		theme(plot.title = element_text(size=9), axis.title = element_text(size=8),axis.text = element_text(size=8))

}

combined <- ggm5[[1]]+ggm5[[2]]+ggm5[[3]]+ggm7[[1]]+ggm7[[2]]+ggm7[[3]] & theme(legend.position="bottom")


pdf(here('figures','si','residuals_lisa_m5_m7_9k_7k.pdf'),width=7.25,height=4.2)
combined+plot_layout(guides='collect',ncol=3)
dev.off()
