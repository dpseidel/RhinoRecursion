### Code for Data parsing & Prep associated with Seidel et al. 2019 
### "Mesoscale movement and recursion behaviors of Namibian black rhinos"
### (submitted to Movement Ecology, BioMove Special Edition)
### Author: Dana Paige Seidel

library(tidyverse)
library(sf)
library(lubridate)
library(magrittr)
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")

# Metadata regarding file index and location of MODIS imagery
# MODIS imagery is freely available from NASA.gov
filemetadata <- read_csv("Data/MODIS_metadata.csv")

# downloaded from http://ecoregions2017.appspot.com/, Licensed under CC-BY 4.0
ecoregions <- read_sf("Data/Ecoregions2017.shp")

# functional water data - available upon reasonable request
# not public in order to protect endangered species
water_utm <- st_read("Data/functional water.shp", crs = 4326) %>%
  st_transform(32733)

# rhino data - available upon reasonable request
# not public in order to protect endangered species
rhinos_pre2017 <- read_sf("Data/clean_rhino_points.shp")
rhinos_1718 <- read_csv("Data/clean_2017-18.csv") %>%
  mutate(
    id = (Tag %>% str_remove(" ") %>% str_remove("IR-")),
    local_date = strftime(with_tz(timestamp, tzone = "Africa/Windhoek"), format = "%m/%d/%y"),
    local_time = strftime(with_tz(timestamp, tzone = "Africa/Windhoek"), format = "%T")
  ) %>%
  st_as_sf(., coords = c("Longitude", "Latitude"), na.fail = FALSE, crs = 4326) %>%
  st_transform(32733) %>%
  mutate(
    x = st_coordinates(.)[, 1],
    y = st_coordinates(.)[, 2]
  ) %>%
  dplyr::select(names(rhinos_pre2017))

# combine datasets
all_rhinos <- rbind(rhinos_pre2017, rhinos_1718) %>%
  mutate(
    date = mdy_hms(paste(local_date, local_time)),
    local_date = mdy(local_date)
  ) %>% # filter 4 white rhinos + 679 which only has 1 pt from future analyses
  filter(!(id %in% c("SAT469", "SAT676", "SAT677", "SAT678", "SAT679"))) 


# Custom function to build a extent/bounding-box polygon from a simple features dataset
sf_polyfrombbox <- function(sf, buffer = 0) {
  if (buffer == 0) {
    bbox <- st_bbox(sf) %>% as.numeric()
  } else {
    bbox0 <- st_bbox(sf) %>% as.numeric()
    bbox <- bbox0 + c(-buffer, -buffer, buffer, buffer)
  }
  
  coords <- list(matrix(c(
    bbox[1], bbox[2],
    bbox[3], bbox[2],
    bbox[3], bbox[4],
    bbox[1], bbox[4],
    bbox[1], bbox[2]
  ),
  ncol = 2, byrow = TRUE
  ))
  st_polygon(coords) %>% st_sfc(crs = st_crs(sf))
}

# transform rhinos to latlong, way faster than projecting the ecoregions shpfile
lat_rhinos <- all_rhinos %>% st_transform(4326)
# Crop the ecoregions shapefile with a small buffer
namib_eco <-  st_intersection(ecoregions,
                              sf_polyfrombbox(st_bbox(lat_rhinos),
                                              buffer = .25)) # ~50km buffer