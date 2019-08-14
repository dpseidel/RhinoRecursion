### Code for 16 day recursion analysis associated with Seidel et al. 2019
### "Mesoscale movement and recursion behaviors of Namibian black rhinos"
### (submitted to Movement Ecology, BioMove Special Edition)
### Author: Dana Paige Seidel

# tlocoh packages are available through R-Forge
# install.packages("tlocoh", repos = "http://R-Forge.R-project.org")
# install.packages("tlocoh.dev", repos = "http://R-Forge.R-project.org")

library(furrr) # for parallel processing
library(tlocoh)
library(tlocoh.dev)
library(velox)
library(nlme)
library(lme4)

source("DataParsing_Script.R") # defines filemetadata, all_rhinos 

########## Data Prep ###############
filemetadata %>% pull(jday) %>% unique() %>% as.numeric() %>% sort()
Interval_Starts <- sort(unique(as.numeric(filemetadata$jday)))

# getting coverage estimates -- on 16 day intervals, require at least 1 fix per day.
filter_df <- all_rhinos %>%
  st_set_geometry(NULL) %>%
  mutate(
    jday = yday(date),
    interval_start = as.numeric(as.character(cut(jday,
      breaks = c(Interval_Starts, 367), # want that last interval to include values through leap years.
      right = F, labels = Interval_Starts
    ))),
    month = lubridate::month(date),
    day = lubridate::day(date),
    year = lubridate::year(date)
  ) %>%
  select(id, interval_start, year, month, day) %>%
  unique() %>%
  group_by(id, interval_start, year) %>%
  tally() %>%
  ungroup() %>%
  mutate(comp90 = case_when(
    interval_start == "353" ~ n / 13,
    year == 2012 & interval_start == "353" ~ n / 14,
    TRUE ~ n / 16
  )) %>%
  filter(comp90 > .90) %>%
  select(ID = id, int = interval_start, Y = year)


rhino_ints <- all_rhinos %>% mutate(
  jday = yday(date),
  interval_start = as.numeric(as.character(cut(jday,
    breaks = c(Interval_Starts, 367), # want that last interval to include values through leap years.
    right = F, labels = Interval_Starts
  )))
)

rhino_trim <- pmap(filter_df, function(ID, int, Y) {
  filter(rhino_ints, id == ID, interval_start == int, year(date) == Y)
})

names(rhino_trim) <- glue::glue_data(filter_df, "{ID}-{Y}-{int}")

dfs <- map(rhino_trim, ~ st_set_geometry(.x, NULL))

extract_NDVI <- function(polys, file) {
  vlx <- velox(file)
  polys <- st_cast(polys, "MULTIPOLYGON") # just to prevent multi-type frames
  polys %>%
    st_set_geometry(NULL) %>%
    mutate(
      meanNDVI =
        vlx$extract(polys, fun = function(x) mean(x, na.rm = TRUE)),
      sdNDVI =
        vlx$extract(polys, fun = function(x) sd(x, na.rm = TRUE)),
      maxNDVI =
        vlx$extract(polys, fun = function(x) max(x, na.rm = TRUE)),
      minNDVI =
        vlx$extract(polys, fun = function(x) min(x, na.rm = TRUE)),
      medianNDVI =
        vlx$extract(polys, fun = function(x) median(x, na.rm = TRUE))
    )
}


########  16 day Recursion Analysis ########

tumaps <- map(dfs, safely(function(df) {
  df <- na.omit(df)
  k <- round(sqrt(nrow(df)))
  lxy <- xyt.lxy(
    xy = matrix(c(df$x, df$y), ncol = 2),
    dt = df$date,
    id = df$id,
    proj4string = CRS("+init=epsg:32733")
  )
  lxy.tumap(lxy, ivg = 12 * 3600, grid = "square", cellsize = 1000)
}))

# # 12 errors: some just too gappy, so we are going to drop them
# compact(map(tumaps, ~.$error))  %>% length()
# compact(map(tumaps, ~.$error)) %>% names()
# # [1] "SAT237-2012-65"   "SAT239-2012-193"  "SAT2590-2018-177" "SAT278-2012-273"  "SAT279-2013-49"   "SAT428-2012-241" 
# # [7] "SAT451-2012-353"  "SAT454-2013-65"   "SAT643-2013-97"   "SAT683-2013-177"  "SAT821-2014-1"    "SAT821-2013-305" 
# drop gappy rhinos
tumaps_sf <- compact(map(tumaps, ~ .$result)) %>%
  map(., ~ .x[[1]] %>% st_as_sf())

tumaps_sf <- pmap(
  list(x = tumaps_sf, y = names(tumaps_sf)),
  function(x, y) {
    split <- str_split(y, "-")[[1]]
    mutate(x, id = split[1], year = split[2], date_id = split[3])
  }
)

polymeta <- names(tumaps_sf) %>%
  str_split("-", simplify = T) %>%
  as_data_frame() %>%
  rename(id = V1, year = V2, date_id = V3) %>%
  mutate(
    file = paste0( # match processed MODIS file names
      "~/Dropbox/Processed/", year,
      formatC(as.numeric(date_id), width = 3, flag = "0"), ".tif"
    ),
    list_name = names(tumaps_sf)
  ) %>%
  arrange(file)

# order to match
ordered_files <- polymeta[match(names(tumaps_sf), polymeta$list_name), ]
# check
# names(tumaps_sf) == ordered_files$list_name # yay! the same!

plan(multiprocess)
biweekly <- future_pmap_dfr(
  list(polys = tumaps_sf, file = ordered_files$file),
  extract_NDVI
)

### Summary Statistics
biweekly %>% select(id, date_id) %>% distinct() %>% nrow()
biweekly %>% select(id) %>% distinct() %>% nrow()
biweekly %>% filter(nsv.43200 > 0) %>% nrow()
biweekly %>%
  filter(nsv.43200 > 0) %>%
  summarise(
    mean_nsv = mean(nsv.43200), sd_nsv = sd(nsv.43200),
    mean_mlnv = mean(mnlv.43200), sd_mnlv = sd(mnlv.43200),
    mean_ndvi = mean(meanNDVI * .0001, na.rm = T), 
    sd_NDVI = sd(sdNDVI * .0001, na.rm = T)
  )

########  16 day Home Range constuction & NDVI Analysis ########
proj4 = "+proj=utm +zone=33 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

construct_UDs <- function(df) {
  ## k-locoh Code
  dropNA <- na.omit(df)
  k <- round(sqrt(nrow(dropNA)))
  lxy <- tlocoh::xyt.lxy(
    xy = matrix(c(dropNA$x, dropNA$y), ncol = 2),
    dt = dropNA$date, id = dropNA$id, tau.diff.max = 0, # disable filter to avoid errors
    proj4string = sp::CRS(proj4),
    status = F
  )
  lxy <- tlocoh::lxy.nn.add(lxy, s = 0, k = k, status = F)
  lhs <- tlocoh::lxy.lhs(lxy,
                         k = k, s = 0, iso.levels = c(0.90),
                         iso.add = T, status = F
  )
  klocoh.isopolys <- tlocoh::isopleths(lhs)[[1]][1:3] %>% st_as_sf()
  
  area_df <- tibble(
    id = df$id[1],
    date_id = paste0(year(df$date[1]), formatC(df$interval_start[1], width = 3, flag = "0")),
    iso.level = c(90),
    method = c("klocoh")
  ) %>%
    mutate(
      area = klocoh.isopolys$area,
      edge = klocoh.isopolys$edge.len
    )
  
  return(list(
    "AreaSummary" = area_df,
    "klocoh" = klocoh.isopolys
  ))
}

plan(multiprocess)
polys <- future_map(dfs, safely(construct_UDs), .progress = T)


join_constructs <- function(x){
  x$result$klocoh %>% 
    select(name = iso.level) %>% 
    mutate(method = "klocoh",
           id = x$result$AreaSummary$id, 
           date_id = x$result$AreaSummary$date_id, 
           iso.level = x$result$AreaSummary$iso.level,
           area = x$result$AreaSummary$area,
           edge = x$result$AreaSummary$edge) %>% 
    select(id, date_id, method, iso.level, area, edge)
}

polys_combined <- map(polys, join_constructs) %>% reduce(rbind)

polys_by_interval <- polys_combined %>% split(.$date_id)

polymeta <- names(polys) %>%
  str_split("-", simplify = T) %>%
  as_tibble() %>%
  rename(id = V1, year = V2, date_id = V3) %>%
  mutate(
    interval = paste0( year,
    formatC(as.numeric(date_id), width = 3, flag = "0")),
    file = paste0(
      "~/Dropbox/Processed/", year,
      formatC(as.numeric(date_id), width = 3, flag = "0"), ".tif"
    ),
    list_name = names(polys)
  ) %>%
  arrange(interval)

unique(polymeta$interval) == names(polys_by_interval)

extracted_UDs <- pmap_df(list(polys = polys_by_interval, 
                             file = unique(polymeta$file)), 
                        extract_NDVI) 

scaledk <- extracted_UDs %>% 
  mutate(iso.level = as.factor(iso.level), 
         date_id = as.character(date_id),
         meanNDVI = meanNDVI *.0001,
         logarea = log(area))%>% 
  rename(met = method)

scaled_k90 <- scaledk %>% filter(met == "klocoh", iso.level == 90)

GLM <- gls(logarea ~ meanNDVI, data = scaled_k90, method = "ML", na.action = na.omit)
lmm <- lme(logarea ~ meanNDVI, data = scaled_k90,
            random = ~1|id, method = "ML",  na.action = na.omit)

anova(GLM, lmm)
summary(lmm)

