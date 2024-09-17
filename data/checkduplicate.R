library(here)
library(rcarbon)
library(geosphere)
library(dplyr)

#Read data ----
d <- read.csv(here('data','c14data.csv')) 
d2  <- select(d,Sitename,Site_code,Latitude,Longitude) |> unique()

# Check for Duplicates in Sitename
any(table(d2$Sitename)>1) #Passed
# Check for Duplicates in Site_code
any(table(d2$Site_code)>1) #Failed

# Identify sitename/sitecode combinations with issues:
d2$ToCheck <- FALSE
d2$ToCheck[d2$Site_code %in% names(which(table(d2$Site_code)>1))] <- TRUE
d2 <- d2[order(d2$Site_code),]
write.csv(d2,'tocheck.csv',row.names=F)
