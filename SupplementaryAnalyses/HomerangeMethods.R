##### LoCoh v. AKDE #####

### To Compare the AKDE vs. k-LoCoH home range 

# replacing the construct_UDs function in 16dayRecursion_HomeRange_Script
# we can extract both types of home range at the same time. 
proj4 = "+proj=utm +zone=33 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

construct_UDs <- function(df) {
  ## CTMM code
  proj4 <- "+proj=utm +zone=33 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
  telm <- stmove:::create_telemetry(df, proj4)
  m.ouf <- ctmm::ctmm.guess(telm, interactive = FALSE) # automated model guess
  M.OUF <- ctmm::ctmm.fit(telm, m.ouf) # this can take awhile...
  UD <- ctmm::akde(telm, M.OUF, weights = F)
  akde.isopolys <- ctmm::SpatialPolygonsDataFrame.UD(UD, level.UD = c(.9)) %>% st_as_sf()
  
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
    iso.level = 90,
    method = c("akde", "klocoh")
  ) %>%
    mutate(
      area = ifelse(method == "klocoh",
                    klocoh.isopolys$area,
                    st_area(akde.isopolys[str_detect(akde.isopolys$name, "ML"), ])
      ),
      edge = ifelse(method == "klocoh", klocoh.isopolys$edge.len,
                    st_length(st_cast(akde.isopolys[str_detect(akde.isopolys$name, "ML"), ], "MULTILINESTRING"))
      )
    )
  
  return(list(
    "AreaSummary" = area_df,
    "AKDE" = akde.isopolys,
    "klocoh" = klocoh.isopolys
  ))
}

plan(multiprocess)
polys <- future_map(dfs, safely(construct_UDs), .progress = T)

join_constructs <- function(x){
  akde <- x$result$AKDE[2,] %>% # ML 90%
    dplyr::mutate(name = as.character(name),
                  method = "akde")
  klocoh <- x$result$klocoh %>% 
    select(name = iso.level) %>% 
    mutate(method = "klocoh")
  
  rbind(akde, klocoh) %>% 
    dplyr::mutate(id = x$result$AreaSummary$id, 
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

extract_NDVI <- function(polys, file){
  vlx <- velox(file)
  polys <- st_cast(polys, "MULTIPOLYGON") # just to prevent multi-type frames
  polys %>% 
    st_set_geometry(NULL) %>% 
    mutate(meanNDVI = 
             vlx$extract(polys, fun = function(x) mean(x, na.rm = TRUE)),
           sdNDVI =
             vlx$extract(polys, fun = function(x) sd(x, na.rm = TRUE)),
           maxNDVI = 
             vlx$extract(polys, fun = function(x) max(x, na.rm = TRUE)),
           minNDVI=
             vlx$extract(polys, fun = function(x) min(x, na.rm = TRUE)),
           medianNDVI = 
             vlx$extract(polys, fun = function(x) median(x, na.rm = TRUE)))
}

extracted_UDs <- pmap_df(list(polys = polys_by_interval, 
                              file = unique(polymeta$file)), 
                         extract_NDVI) 

scaled <- extracted_UDs %>% 
  mutate(iso.level = as.factor(iso.level), 
         date_id = as.character(date_id),
         meanNDVI = meanNDVI *.0001,
         logarea = log(area))%>% 
  rename(met = method)


scaled_k90 <- scaled %>% filter(met == "klocoh", iso.level == 90)
scaled_a90 <- scaled %>% filter(met == "akde", iso.level == 90)


lmm <- lme(logarea ~ meanNDVI, data = scaled_k90,
           random = ~1|id, method = "ML",  na.action = na.omit)

summary(lmm)

lmm_a90 <- lme(logarea ~ meanNDVI, data = scaled_a90,
           random = ~1|id, method = "ML",  na.action = na.omit)
summary(lmm_a90)
