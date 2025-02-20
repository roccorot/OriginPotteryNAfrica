# Load required libraries
library(ggplot2)
library(sf)
library(readr)
library(readxl)
library(rnaturalearth)
library(dplyr)
library(elevatr)
library(raster)
library(cowplot)
library(ggrepel)
library(marmap)
library(here)

# Load data
sites <- read_csv("c14data.csv")
key_sites <- read_csv("Key sites.csv")

# Remove duplicate sites based on latitude and longitude
sites_unique <- sites %>% distinct(Latitude, Longitude, .keep_all = TRUE)
key_sites_unique <- key_sites %>% distinct(Latitude, Longitude, .keep_all = TRUE)

# Order key sites alphabetically for legend numbering
key_sites_unique <- key_sites_unique %>% arrange(Sitename)
key_sites_unique$Label <- seq_along(key_sites_unique$Sitename)

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
  annotate("point", x = -20, y = 13, shape = 20, size = 2, color = "black") +
  annotate("text", x = -19, y = 13, label = "Archaeological Sites", hjust = 0, size = 2, fontface = "bold") +
  annotate("point", x = -20, y = 12, shape = 24, size = 2, fill = "red", color = "black") +
  annotate("text", x = -19, y = 12, label = "Key Sites", hjust = 0, size = 2, fontface = "bold") +
  annotate("text", x = -20, y = 40, label = legend_text, hjust = 0, vjust = 1, 
           size = 2, fontface = "bold", color = "black")

# Export to PDF with proper size
  ggsave("Figure1.pdf", plot = fig1, width = 7.25, height = 4.2, units = "in", dpi = 300)
