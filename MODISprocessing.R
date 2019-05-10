### Code for MODIS Imagery Processing associated with Seidel et al. 2019
### "Mesoscale movement and recursion behaviors of Namibian black rhinos"
### (submitted to Movement Ecology, BioMove Special Edition)
### Author: Dana Paige Seidel

library(gdalUtils)
library(tidyverse)
library(lubridate)

# MODIS imagery downloaded from NASA earth explorer
# Imagery is downloaded as HDF files and must be converted to geotiffs, tiles 
# mosaicked, and the images projected and cropped by extent of rhinos

hdf_path <- "MODIS/Raw"
hdf_files <- list.files(hdf_path, full.names = T, pattern = ".hdf$")

# some gaps due to cloud cover... need to pull those files too. 9 files
cloudfiles <- list.files(paste0(hdf_path, "/cloudcover_g10"), full.names = T, pattern = ".hdf$")

julian2date <- function(year, jday) {
  x <- as.Date(paste(year, "01-01", sep = "-"))
  yday(x) <- as.numeric(jday)
  return(x)
}

filemetadata <- tibble(
  path = c(hdf_files, cloudfiles),
  group = str_extract(path, pattern = "v\\d{2}"),
  index = str_extract(path, pattern = "A\\d{7}"),
  year = str_extract(index, "\\d{4}"),
  jday = str_extract(index, "\\d{3}$"),
  start_date = julian2date(year, jday),
  end_date = start_date + 15,
  notes = ifelse(str_detect(path, "cloud"), "cloud cover greater than 10% in raw image.", NA)
) %>%
  arrange(start_date)

Interval_Starts <- sort(unique(as.numeric(filemetadata$jday)))

# mosaic north and south images, project, and crop
prepMODIS <- function(files, temp_dir = tempdir()) {
  batch_gdal_translate(files,
    sd_index = 1,
    outdir = temp_dir,
    outsuffix = ".tif"
  )

  file_groups <- tibble(
    file = list.files(temp_dir, full.names = T),
    date_id = str_extract(file, "\\d{7}")
  ) %>%
    split(.$date_id)

  walk(file_groups, function(files) {
    gdalUtils::mosaic_rasters(files$file, dst_dataset = paste0(temp_dir, files$date_id[1], ".tif"))
  })

  walk(names(file_groups), ~ gdalwarp(paste0(temp_dir, .x, ".tif"),
    paste0("MODIS/Processed/", .x, ".tif"), # requires absolute path
    s_srs = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs",
    t_srs = "+proj=utm +zone=33 +south +datum=WGS84 +units=m +no_defs",
    te = as.numeric(st_bbox(cropbox))
  ))

  file.remove(list.files(temp_dir))
}

prepMODIS(c(hdf_files, cloudfiles), temp_dir = "MODIStemp/")
