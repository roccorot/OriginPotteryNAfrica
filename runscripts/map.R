# Load required libraries
library(ggplot2)
library(sf)
library(rnaturalearth)
library(dplyr)
library(elevatr)
library(raster)
library(cowplot)
library(ggrepel)
library(marmap)
library(here)

# Load data
sites <- read.csv(here('data','c14data.csv'))
# Define key sites
key_site_codes <- c('ABUS10','AMK','E798','Bu44','E06-1','ONJ_RDM','SRB2','TGL','TAK','TMY3','TMT','TTE','UAF','UTB','WEA83/33')
key_site_names <- c('Adrar Bous 10', 'Amekni', 'Bir Kiseiba E79-8', 'Busharia I', 'Nabta Playa E06-1', 'Ounjougou', 'Sarourab II', 'Tagalagal', 'Takarkori', 'Tamaya Mellet TMY3', 'Temet', 'Tin Torha East', 'Uan Afuda', 'Uan Tabu', 'Wadi el Akhdar 83/33')

# Remove duplicate sites based on latitude and longitude
# sites_unique <- sites %>% distinct(Latitude, Longitude, .keep_all = TRUE)
sites_unique <- unique(sites[,c('Site_code','Latitude','Longitude')])
key_sites_unique <- subset(sites_unique,Site_code%in%key_site_codes)
key_sites_unique$Sitename <- key_site_names
key_sites_unique$Label <- 1:nrow(key_sites_unique)

# Convert data frames to sf objects
sites_sf <- st_as_sf(sites_unique, coords = c("Longitude", "Latitude"), crs = 4326)
key_sites_sf <- st_as_sf(key_sites_unique, coords = c("Longitude", "Latitude"), crs = 4326)

# Get North Africa map
north_africa <- ne_countries(continent = "Africa", scale = "medium", returnclass = "sf")

# Get elevation data (Restricting to North Africa bounding box)
elevation <- get_elev_raster(locations = north_africa, z = 3, clip = "locations", src = 'aws')

# Convert elevation data to dataframe
dem <- as.data.frame(elevation, xy = TRUE) |> na.omit()
colnames(dem)[3] <- "elevation"
dem$elevation[dem$elevation < 0] <- 0  # Set sea level elevations to 0


# Generate the map
fig1 <- ggplot() +
  geom_tile(data = dem, aes(x = x, y = y, fill = elevation)) +
  scale_fill_etopo(guide="none") +
  geom_sf(data = sites_sf, size = 1.5, col = 'black', pch = 20) +
  geom_sf(data = key_sites_sf, color = "black", shape = 24, fill = "red", size = 3) +
  geom_label_repel(data = key_sites_unique, aes(x = Longitude, y = Latitude, label = Label), 
                   size = 3, box.padding = 0.5, point.padding = 0.3, segment.color = "black")+
  coord_sf(xlim = c(-20, 40), ylim = c(10, 40)) +  
  theme_minimal() +
  labs(title = "", x= "", y="")


# Create legend text
legend_data <- key_sites_unique %>% dplyr::select(Label, Sitename)
legend_text <- paste(legend_data$Label, legend_data$Sitename, sep = " - ")
legend_text <- paste(legend_text, collapse = "\n")

# Add legend 

  fig1 <- fig1 +
  annotate("point", x = -20, y = 13, shape = 20, size = 2.5, color = "black") +
  annotate("text", x = -19, y = 13, label = "Archaeological Sites", hjust = 0, size = 2.5, fontface = "bold") +
  annotate("point", x = -20, y = 12, shape = 24, size = 2.5, fill = "red", color = "black") +
  annotate("text", x = -19, y = 12, label = "Key Sites", hjust = 0, size = 2.5, fontface = "bold") +
  annotate("text", x = -20, y = 40, label = legend_text, hjust = 0, vjust = 1, 
           size = 2.5, fontface = "bold", color = "black")

# Export to PDF with proper size
  ggsave(here('figures','main','map.pdf'), plot = fig1, width = 7.25, height = 4.2, units = "in", dpi = 300)
