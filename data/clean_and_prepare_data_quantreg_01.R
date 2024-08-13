#Load relevant libraries ----
library(here)
library(rcarbon)
library(geosphere)
library(dplyr)

#Read data ----
d <- read.csv(here('data','c14data.csv'))

# Calibrate, Bin, and Thin ----
dcal <- calibrate(d$Yrs.BP,d$STD) #calibrate
d$medcal <- medCal(dcal) #store median calibrated date
d$theta <- d$medcal #define initialisation of theta
d$SiteID  <- as.numeric(as.factor(d$Sitename)) #assign unique SiteID as integer
d$bin <- binPrep(sites=d$SiteID,ages=dcal,h=100)  #binning with 100yrs interval
bin.table <- aggregate(pottery~bin,max,data=d) |> setNames(c('bin','pottery_bin')) #aggegate pottery presence/absence by bin
thinned <- thinDates(ages=d$Yrs.BP,errors=d$STD,bins=d$bin,size=1,thresh=1,seed=123,method='splitsample') #random thinning
d <- d[thinned,]
d  <- left_join(d,bin.table)
d  <- subset(d,pottery_bin==1) #Subset to context/sites with presence only
# Redate and re-assign IDs
d$DateID <- 1:nrow(d)
d$SiteID  <- as.numeric(as.factor(d$Sitename)) #assign unique SiteID as integer
dcal <- calibrate(d$Yrs.BP,d$STD) #calibrate
theta <- d$theta
# Handle Stratigraphic Relationship
site_subset <- select(d,SiteID,STRATIGRAPHY) |> unique() |> subset(STRATIGRAPHY!='')
relationships <- numeric()

for (i in 1:nrow(site_subset))
{
	tmp <- site_subset$STRATIGRAPHY[i]
	extracted_strings <- regmatches(tmp, gregexpr("\\[(.*?)\\]", tmp, perl=TRUE))[[1]] |> gsub(pattern="\\[|\\]", replacement="",x=_)
	d2  <- subset(d,SiteID==site_subset$SiteID[i])
	for (j in (1:length(extracted_strings)-1))
	{
		left <- d2$DateID[d2$Layer.SU.UNIT.LOCUS==extracted_strings[j]]
		right <- d2$DateID[d2$Layer.SU.UNIT.LOCUS==extracted_strings[j+1]]
		left <- left[which(!is.na(left))]
		right <- right[which(!is.na(right))]
		if (any(is.na(left))){print(i)}
		if(length(left)==1 & length(right)==1)
		{
			relationships <- c(relationships,paste0('theta[',left,']>theta[',right,']'))
		}
		if (length(left)>=1 & length(right>=1) & (length(left)>1 | length(right)>1))
		{
			for (k1 in 1:length(left))
			{
				for (k2 in 1:length(right))
				{
					relationships <- c(relationships,paste0('theta[',left[k1],']>theta[',right[k2],']'))
				}
			}
		}	
	}

}

check  <- unlist(sapply(relationships,FUN=function(x){eval(parse(text=x))}))

 while(!all(check))
 {
	for (i in 1:length(relationships))
	{
		indices <- regmatches(relationships[i], gregexpr("\\[(.*?)\\]", relationships[i], perl=TRUE))[[1]] |> gsub(pattern="\\[|\\]", replacement="",x=_) |> as.numeric()
		if (theta[indices[1]] <= theta[indices[2]])
		{
			pool1 <- dcal$grids[[indices[1]]]
			pool2 <- dcal$grids[[indices[2]]]
			while(theta[indices[1]] <= theta[indices[2]])
			{
				theta[indices[1]] <- sample(pool1$calBP,size=1,prob=pool1$PrDens)
				theta[indices[2]] <- sample(pool2$calBP,size=1,prob=pool2$PrDens)
			}
		}
	}
	check  <- unlist(sapply(relationships,FUN=function(x){eval(parse(text=x))}))
}

d$theta <- theta #re-assign initial theta values
d$theta.scaled <- scale(d$theta) |> as.numeric()
theta.constraints <- paste(relationships,collapse=' & ')


# Compute Distances from Putative Origins ----
# Compute distance from origins
d$Latitude <- round(d$Latitude,2)
d$Longitude <- round(d$Longitude,2)
e798 <- d[which(d$Site_code=='E798')[1],c('Longitude','Latitude')]
abus10 <- d[which(d$Site_code=='ABUS10')[1],c('Longitude','Latitude')]
onj_rdm <- d[which(d$Site_code=='ONJ_RDM')[1],c('Longitude','Latitude')]

d$dist.e <- distVincentyEllipsoid(e798,d[,c('Longitude','Latitude')]) |> round()
d$dist.a <- distVincentyEllipsoid(abus10,d[,c('Longitude','Latitude')]) |> round()
d$dist.o <- distVincentyEllipsoid(onj_rdm,d[,c('Longitude','Latitude')]) |> round()

# Record stratigraphic relations 
constraint.text  <- paste0('chronological.constraint ~ dconstraint(',theta.constraints,')')

save(d,constraint.text,file=here('data','input_quantreg.RData'))
