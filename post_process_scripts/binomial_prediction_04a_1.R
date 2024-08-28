# Load Libraries ----
library(sf)
library(patchwork)
library(terra)
library(rcarbon)
library(raster)
library(dplyr)
library(here)
library(ggplot2)
library(rnaturalearth)

# Load Observed Data and Results ----
load(here('data','input_binom.RData'))
load(here('results_images','binom','results_binom_m1.RData'))
load(here('results_images','binom','results_binom_m2.RData'))
load(here('results_images','binom','results_binom_m3.RData'))
load(here('results_images','binom','results_binom_m4.RData'))
load(here('results_images','binom','results_binom_m5.RData'))
load(here('results_images','binom','results_binom_m6.RData'))
load(here('results_images','binom','results_binom_m7.RData'))

# Read window of analyses ----
win <- st_read(here('data','shp','Research_area.shp')) |> st_geometry()
world <- ne_countries(scale = "medium", returnclass = "sf")
background <- st_geometry(world)

# Rasterise and compute distance from origin(s) ----
template  <-  rast(vect(win),res=0.5)
win <- rasterize(vect(win), template)

# coordinates in long - lat
e798 <- d[which(d$Site_code=='E798')[1],c('Longitude','Latitude')] |> as.numeric()
abus10 <- d[which(d$Site_code=='ABUS10')[1],c('Longitude','Latitude')] |> as.numeric()
onj_rdm <- d[which(d$Site_code=='ONJ_RDM')[1],c('Longitude','Latitude')] |> as.numeric()

d.e <- pointDistance(matrix(e798,ncol=2),matrix(crds(win,na.rm=F),ncol=2),longlat=TRUE)
d.e[which(is.na(values(win)))] <- NA
d.a <- pointDistance(matrix(abus10,ncol=2),matrix(crds(win,na.rm=F),ncol=2),longlat=TRUE)
d.a[which(is.na(values(win)))] <- NA
d.o <- pointDistance(matrix(onj_rdm,ncol=2),matrix(crds(win,na.rm=F),ncol=2),longlat=TRUE)
d.o[which(is.na(values(win)))] <- NA

win.o <- win.a <- win.e <- win
values(win.o) <- d.o
values(win.a) <- d.a
values(win.e) <- d.e

raster.e <- as.data.frame(win.e,xy=T)
raster.a <- as.data.frame(win.a,xy=T)
raster.o <- as.data.frame(win.o,xy=T)
raster_predictions <- data.frame(x=raster.e$x,y=raster.e$y,dist.e=raster.e$layer,dist.a=raster.a$layer,dist.o=raster.o$layer)

mean.dists <- mean(c(d$dist.e,d$dist.a,d$dist.o))
sd.dists <- sd(c(d$dist.e,d$dist.a,d$dist.o))
raster_predictions$dist.e.scaled <- (raster_predictions$dist.e - mean.dists)/sd.dists
raster_predictions$dist.a.scaled <- (raster_predictions$dist.a - mean.dists)/sd.dists
raster_predictions$dist.o.scaled <- (raster_predictions$dist.o - mean.dists)/sd.dists

# Fill predictions ----
time.slices  <- seq(12000,7000,-1000)
time.slices.scaled <- (time.slices - mean(d$theta))/sd(d$theta)
nsim <- 1000
raster_predictions$m1 <- raster_predictions$m2 <- raster_predictions$m3 <- raster_predictions$m4 <- raster_predictions$m5 <- raster_predictions$m6 <- raster_predictions$m7 <- NA

results <- vector('list',length=length(time.slices))
sites <- vector('list',length=length(time.slices))
d.calibrated  <- calibrate(d$Yrs.BP,d$STD)
d$pottery_bin_lab <- ifelse(d$pottery_bin==1,'Present','Absent')
d$pottery_bin_lab <- factor(d$pottery_bin_lab)
sites[[1]] <- st_as_sf(d[which.CalDates(d.calibrated,BP < (12000+500) & BP > (12000 - 500),p=0.5),],coords=c('Longitude','Latitude'),crs=4326)
sites[[2]] <- st_as_sf(d[which.CalDates(d.calibrated,BP < (11000+500) & BP > (11000 - 500),p=0.5),],coords=c('Longitude','Latitude'),crs=4326)
sites[[3]] <- st_as_sf(d[which.CalDates(d.calibrated,BP < (10000+500) & BP > (10000 - 500),p=0.5),],coords=c('Longitude','Latitude'),crs=4326)
sites[[4]] <- st_as_sf(d[which.CalDates(d.calibrated,BP < (9000+500) & BP > (9000 - 500),p=0.5),],coords=c('Longitude','Latitude'),crs=4326)
sites[[5]] <- st_as_sf(d[which.CalDates(d.calibrated,BP < (8000+500) & BP > (8000 - 500),p=0.5),],coords=c('Longitude','Latitude'),crs=4326)
sites[[6]] <- st_as_sf(d[which.CalDates(d.calibrated,BP < (7000+500) & BP > (7000 - 500),p=0.5),],coords=c('Longitude','Latitude'),crs=4326)

set.seed(123)
sim.index <- sample(1:nrow(m1.res$post),size=1000)

for(i in 1:length(time.slices))
{
	results[[i]] <- raster_predictions 
	time.x  <- time.slices.scaled[i]

	#m1
	post.m1 <- m1.res$post[sim.index,]
	results[[i]]$m1 <- sapply(1:nsim,function(x,post,t,d){return(1/(1+exp(-(post[x,1] - post[x,2]*t - post[x,3]*d))))},post=post.m1,t=time.x,d=results[[i]]$dist.e.scaled) |> apply(1,mean)
	#m2
	post.m2 <- m2.res$post[sim.index,]
	results[[i]]$m2 <- sapply(1:nsim,function(x,post,t,d){return(1/(1+exp(-(post[x,1] - post[x,2]*t - post[x,3]*d))))},post=post.m2,t=time.x,d=results[[i]]$dist.a.scaled) |> apply(1,mean)
	#m3
	post.m3 <- m3.res$post[sim.index,]
	results[[i]]$m3 <- sapply(1:nsim,function(x,post,t,d){return(1/(1+exp(-(post[x,1] - post[x,2]*t - post[x,3]*d))))},post=post.m2,t=time.x,d=results[[i]]$dist.o.scaled) |> apply(1,mean)
	#m4  
	post.m4 <- m4.res$post[sim.index,]
	results[[i]]$m4 <- sapply(1:nsim,function(x,post,t,d1,d2){return(pmax(c(1/(1+exp(-(post[x,1] - post[x,3]*t - post[x,5]*d1)))),c(1/(1+exp(-(post[x,2] - post[x,4]*t - post[x,6]*d2))))))},post=post.m4,t=time.x,d1=results[[i]]$dist.e.scaled,d2=results[[i]]$dist.a.scaled) |> apply(1,mean)
	#m5  
	post.m5 <- m5.res$post[sim.index,]
	results[[i]]$m5 <- sapply(1:nsim,function(x,post,t,d1,d2){return(pmax(c(1/(1+exp(-(post[x,1] - post[x,3]*t - post[x,5]*d1)))),c(1/(1+exp(-(post[x,2] - post[x,4]*t - post[x,6]*d2))))))},post=post.m5,t=time.x,d1=results[[i]]$dist.e.scaled,d2=results[[i]]$dist.o.scaled) |> apply(1,mean)
	#m6  
	post.m6 <- m6.res$post[sim.index,]
	results[[i]]$m6 <- sapply(1:nsim,function(x,post,t,d1,d2){return(pmax(c(1/(1+exp(-(post[x,1] - post[x,3]*t - post[x,5]*d1)))),c(1/(1+exp(-(post[x,2] - post[x,4]*t - post[x,6]*d2))))))},post=post.m6,t=time.x,d1=results[[i]]$dist.a.scaled,d2=results[[i]]$dist.o.scaled) |> apply(1,mean)
	#m7  
	post.m7 <- m7.res$post[sim.index,]
	results[[i]]$m7 <- sapply(1:nsim,function(x,post,t,d1,d2,d3){return(pmax(c(1/(1+exp(-(post[x,1] - post[x,4]*t - post[x,7]*d1)))),c(1/(1+exp(-(post[x,2] - post[x,5]*t - post[x,8]*d2)))),c(1/(1+exp(-(post[x,3] - post[x,6]*t - post[x,9]*d3))))))},post=post.m7,t=time.x,d1=results[[i]]$dist.e.scaled,d2=results[[i]]$dist.a.scaled,d3=results[[i]]$dist.o.scaled) |> apply(1,mean)
}


ggplots.m1 <- ggplots.m2 <- ggplots.m3 <- ggplots.m4 <- ggplots.m5 <- ggplots.m6 <- ggplots.m7 <- ggplots.m8  <- vector('list',length=length(time.slices))


legend_data <- data.frame(
  presence = factor(c('Present', 'Absent'), levels = c('Present', 'Absent')),
  x = c(0, 0),  # Dummy coordinates
  y = c(0, 0)
)

for(i in 1:length(time.slices))
{

	ggplots.m1[[i]] <- ggplot() +
		geom_sf(data = background) +
		geom_raster(data = results[[i]], aes(x = x, y = y, fill = m1)) +
		geom_sf(data = sites[[i]], aes(shape = factor(pottery_bin)),show.legend=FALSE) +
		geom_point(data = legend_data, aes(x = x, y = y, shape = presence), size = 3, show.legend = TRUE) +  # Add dummy points for legend
		geom_sf(data = background,fill=NA) +
 		scale_fill_viridis_c(limits=c(0,1)) +
		xlim(-17.5,40.5) +
		ylim(7.5,34.5) +
		scale_shape_manual(limits=c('0','1'),values=c('0'=1,'1'=19),labels=c('1'='Present','0'='Absent'),drop=F,guide = guide_legend(override.aes=list(shape=c(1,19)))) +
		labs(title = paste0('BP ',time.slices[i]), x = "Longitude",   y = "Latitude", shape='Pottery',fill='Probability') +
		theme_minimal() +
		theme(plot.title = element_text(size=9), axis.title = element_text(size=8),axis.text = element_text(size=8))

	ggplots.m2[[i]] <- ggplot() +
		geom_sf(data = background) +
		geom_raster(data = results[[i]], aes(x = x, y = y, fill = m2)) +
		geom_sf(data = sites[[i]], aes(shape = factor(pottery_bin)),show.legend=FALSE) +
		geom_point(data = legend_data, aes(x = x, y = y, shape = presence), size = 3, show.legend = TRUE) +  # Add dummy points for legend
		geom_sf(data = background,fill=NA) +
 		scale_fill_viridis_c(limits=c(0,1)) +
		xlim(-17.5,40.5) +
		ylim(7.5,34.5) +
		scale_shape_manual(limits=c('0','1'),values=c('0'=1,'1'=19),labels=c('1'='Present','0'='Absent'),drop=F,guide = guide_legend(override.aes=list(shape=c(1,19)))) +
		labs(title = paste0('BP ',time.slices[i]), x = "Longitude",   y = "Latitude", shape='Pottery',fill='Probability') +
		theme_minimal() +
		theme(plot.title = element_text(size=9), axis.title = element_text(size=8),axis.text = element_text(size=8))

	ggplots.m3[[i]] <- ggplot() +
		geom_sf(data = background) +
		geom_raster(data = results[[i]], aes(x = x, y = y, fill = m3)) +
		geom_sf(data = sites[[i]], aes(shape = factor(pottery_bin)),show.legend=FALSE) +
		geom_point(data = legend_data, aes(x = x, y = y, shape = presence), size = 3, show.legend = TRUE) +  # Add dummy points for legend
		geom_sf(data = background,fill=NA) +
 		scale_fill_viridis_c(limits=c(0,1)) +
		xlim(-17.5,40.5) +
		ylim(7.5,34.5) +
		scale_shape_manual(limits=c('0','1'),values=c('0'=1,'1'=19),labels=c('1'='Present','0'='Absent'),drop=F,guide = guide_legend(override.aes=list(shape=c(1,19)))) +
		labs(title = paste0('BP ',time.slices[i]), x = "Longitude",   y = "Latitude", shape='Pottery',fill='Probability') +
		theme_minimal() +
		theme(plot.title = element_text(size=9), axis.title = element_text(size=8),axis.text = element_text(size=8))


	ggplots.m4[[i]] <- ggplot() +
		geom_sf(data = background) +
		geom_raster(data = results[[i]], aes(x = x, y = y, fill = m3)) +
		geom_sf(data = sites[[i]], aes(shape = factor(pottery_bin)),show.legend=FALSE) +
		geom_point(data = legend_data, aes(x = x, y = y, shape = presence), size = 3, show.legend = TRUE) +  # Add dummy points for legend
		geom_sf(data = background,fill=NA) +
 		scale_fill_viridis_c(limits=c(0,1)) +
		xlim(-17.5,40.5) +
		ylim(7.5,34.5) +
		scale_shape_manual(limits=c('0','1'),values=c('0'=1,'1'=19),labels=c('1'='Present','0'='Absent'),drop=F,guide = guide_legend(override.aes=list(shape=c(1,19)))) +
		labs(title = paste0('BP ',time.slices[i]), x = "Longitude",   y = "Latitude", shape='Pottery',fill='Probability') +
		theme_minimal() +
		theme(plot.title = element_text(size=9), axis.title = element_text(size=8),axis.text = element_text(size=8))

	ggplots.m5[[i]] <- ggplot() +
		geom_sf(data = background) +
		geom_raster(data = results[[i]], aes(x = x, y = y, fill = m3)) +
		geom_sf(data = sites[[i]], aes(shape = factor(pottery_bin)),show.legend=FALSE) +
		geom_point(data = legend_data, aes(x = x, y = y, shape = presence), size = 3, show.legend = TRUE) +  # Add dummy points for legend
		geom_sf(data = background,fill=NA) +
 		scale_fill_viridis_c(limits=c(0,1)) +
		xlim(-17.5,40.5) +
		ylim(7.5,34.5) +
		scale_shape_manual(limits=c('0','1'),values=c('0'=1,'1'=19),labels=c('1'='Present','0'='Absent'),drop=F,guide = guide_legend(override.aes=list(shape=c(1,19)))) +
		labs(title = paste0('BP ',time.slices[i]), x = "Longitude",   y = "Latitude", shape='Pottery',fill='Probability') +
		theme_minimal() +
		theme(plot.title = element_text(size=9), axis.title = element_text(size=8),axis.text = element_text(size=8))

	ggplots.m6[[i]] <- ggplot() +
		geom_sf(data = background) +
		geom_raster(data = results[[i]], aes(x = x, y = y, fill = m3)) +
		geom_sf(data = sites[[i]], aes(shape = factor(pottery_bin)),show.legend=FALSE) +
		geom_point(data = legend_data, aes(x = x, y = y, shape = presence), size = 3, show.legend = TRUE) +  # Add dummy points for legend
		geom_sf(data = background,fill=NA) +
 		scale_fill_viridis_c(limits=c(0,1)) +
		xlim(-17.5,40.5) +
		ylim(7.5,34.5) +
		scale_shape_manual(limits=c('0','1'),values=c('0'=1,'1'=19),labels=c('1'='Present','0'='Absent'),drop=F,guide = guide_legend(override.aes=list(shape=c(1,19)))) +
		labs(title = paste0('BP ',time.slices[i]), x = "Longitude",   y = "Latitude", shape='Pottery',fill='Probability') +
		theme_minimal() +
		theme(plot.title = element_text(size=9), axis.title = element_text(size=8),axis.text = element_text(size=8))

	ggplots.m7[[i]] <- ggplot() +
		geom_sf(data = background) +
		geom_raster(data = results[[i]], aes(x = x, y = y, fill = m7)) +
		geom_sf(data = sites[[i]], aes(shape = factor(pottery_bin)),show.legend=FALSE) +
		geom_point(data = legend_data, aes(x = x, y = y, shape = presence), size = 3, show.legend = TRUE) +  # Add dummy points for legend
		geom_sf(data = background,fill=NA) +
 		scale_fill_viridis_c(limits=c(0,1)) +
		xlim(-17.5,40.5) +
		ylim(7.5,34.5) +
		scale_shape_manual(limits=c('0','1'),values=c('0'=1,'1'=19),labels=c('1'='Present','0'='Absent'),drop=F,guide = guide_legend(override.aes=list(shape=c(1,19)))) +
		labs(title = paste0('BP ',time.slices[i]), x = "Longitude",   y = "Latitude", shape='Pottery',fill='Probability') +
		theme_minimal() +
		theme(plot.title = element_text(size=9), axis.title = element_text(size=8), axis.text = element_text(size=8))
}

combined.m1 <- ggplots.m1[[1]] + ggplots.m1[[2]] + ggplots.m1[[3]] + ggplots.m1[[4]] + ggplots.m1[[5]] + ggplots.m1[[6]] & theme(legend.position = "bottom")
combined.m2 <- ggplots.m2[[1]] + ggplots.m2[[2]] + ggplots.m2[[3]] + ggplots.m2[[4]] + ggplots.m2[[5]] + ggplots.m2[[6]] & theme(legend.position = "bottom")
combined.m3 <- ggplots.m3[[1]] + ggplots.m3[[2]] + ggplots.m3[[3]] + ggplots.m3[[4]] + ggplots.m3[[5]] + ggplots.m3[[6]] & theme(legend.position = "bottom")
combined.m4 <- ggplots.m4[[1]] + ggplots.m4[[2]] + ggplots.m4[[3]] + ggplots.m4[[4]] + ggplots.m4[[5]] + ggplots.m4[[6]] & theme(legend.position = "bottom")
combined.m5 <- ggplots.m5[[1]] + ggplots.m5[[2]] + ggplots.m5[[3]] + ggplots.m5[[4]] + ggplots.m5[[5]] + ggplots.m5[[6]] & theme(legend.position = "bottom")
combined.m6 <- ggplots.m6[[1]] + ggplots.m6[[2]] + ggplots.m6[[3]] + ggplots.m6[[4]] + ggplots.m6[[5]] + ggplots.m6[[6]] & theme(legend.position = "bottom")
combined.m7 <- ggplots.m7[[1]] + ggplots.m7[[2]] + ggplots.m7[[3]] + ggplots.m7[[4]] + ggplots.m7[[5]] + ggplots.m7[[6]] & theme(legend.position = "bottom")

# pdf(here('figures','si','pred_m1b.pdf'),width=10,height=8)
pdf(here('figures','si','pred_m1b.pdf'),width=7.25,height=4.2)
combined.m1 + plot_layout(guides = "collect",ncol=3)
dev.off()

pdf(here('figures','si','pred_m2b.pdf'),width=7.25,height=4.2)
combined.m2 + plot_layout(guides = "collect",ncol=3)
dev.off()

pdf(here('figures','si','pred_m3b.pdf'),width=7.25,height=4.2)
combined.m3 + plot_layout(guides = "collect",ncol=3)
dev.off()

pdf(here('figures','si','pred_m4b.pdf'),width=7.25,height=4.2)
combined.m4 + plot_layout(guides = "collect",ncol=3)
dev.off()

pdf(here('figures','main','pred_m5b.pdf'),width=7.25,height=4.2)
combined.m5 + plot_layout(guides = "collect",ncol=3)
dev.off()

pdf(here('figures','si','pred_m6b.pdf'),width=7.25,height=4.2)
combined.m6 + plot_layout(guides = "collect",ncol=3)
dev.off()

pdf(here('figures','main','pred_m7b.pdf'),width=7.25,height=4.2)
combined.m7 + plot_layout(guides = "collect",ncol=3)
dev.off()
