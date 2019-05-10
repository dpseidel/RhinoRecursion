### Code for 24 hour displacement analysis associated with Seidel et al. 2019 
### "Mesoscale movement and recursion behaviors of Namibian black rhinos"
### (submitted to Movement Ecology, BioMove Special Edition)
### Author: Dana Paige Seidel

source("DataParsing_Script.R")
library(caret)

## Time of Day displacement analysis
ToDfilter <- all_rhinos  %>% st_set_geometry(NULL) %>% 
  mutate(rnd_hour = hour(round_date(date, unit = "1 hour")), 
         ToD = case_when(
           rnd_hour %in% c(23, 0:3) ~ "mid-night",
           rnd_hour %in% 5:9 ~ "dawn",
           rnd_hour %in% 11:15 ~ "mid-day",
           rnd_hour %in% 17:21 ~ "dusk",
           TRUE ~ "other"
         ), 
         ToD_est = case_when(
           rnd_hour %in% c(23, 0:3) ~ 1,
           rnd_hour %in% 5:9 ~ 7 ,
           rnd_hour %in% 11:15 ~ 13,
           rnd_hour %in% 17:21 ~ 19,
           TRUE ~ 0
         ),
         offset_hour = ymd_hms(paste(local_date, local_time)) - 
           ymd_h(paste(local_date, rnd_hour)),
         offset_ToD = 
           ymd_hms(paste(local_date, local_time)) - 
           ymd_h(paste(local_date, ToD_est))) %>% 
  filter(ToD != "other") %>% 
  group_by(id, local_date, ToD) %>% 
  mutate(closest =  abs(as.numeric(offset_ToD)) == min(abs(as.numeric(offset_ToD)))) %>% 
  filter(closest) %>% ungroup() %>% 
  mutate(ToD = readr::parse_factor(ToD, ordered = T, levels = c("mid-night", "dawn", "mid-day", "dusk")))

ToD <- ToDfilter %>% split(.$id) %>% 
  map(., ~.x %>% 
        complete(local_date = seq.Date(from = min(local_date, na.rm = T), to = max(local_date, na.rm = T), 1)) %>% 
        complete(ToD, nesting(local_date)) %>% filter(!is.na(ToD)) %>% 
        arrange(local_date, ToD)) 

lag4 <- ToD %>%
  keep(., function(x) {
    na.omit(unique(x$id)) != "SAT645"
  }) %>% # drop 645, too short, throws error.
  map_df(., function(x) {
    mutate(x,
           dx = c(diff(x), NA),
           dy = c(diff(y), NA),
           dt = sqrt(dx^2 + dy^2),
           dx2 = c(diff(x, lag = 2), NA, NA),
           dy2 = c(diff(y, lag = 2), NA, NA),
           dt2 = sqrt(dx2^2 + dy2^2),
           dx4 = c(diff(x, lag = 4), NA, NA, NA, NA),
           dy4 = c(diff(y, lag = 4), NA, NA, NA, NA),
           dt4 = sqrt(dx4^2 + dy4^2),
           n = n()
    )
  })

# identify individuals with consistent 12 hr fixes 
con_12hr <- map(ToD, function(df) {
  df <- filter(df, ToD %in% c("dawn", "dusk"))
  con <- na.contiguous(df$x)
  df[attr(con, "tsp")[1]:attr(con, "tsp")[2], ]
})

ids <- names(map_dbl(con_12hr, nrow)[map_dbl(con_12hr, nrow) > 60])

chosen9 <- con_12hr[ids] %>%
  map_df(., ~ mutate(.x,
                     dx = c(diff(x), NA),
                     dy = c(diff(y), NA),
                     dt = sqrt(dx^2 + dy^2),
                     dx2 = c(diff(x, lag = 2), NA, NA),
                     dy2 = c(diff(y, lag = 2), NA, NA),
                     dt2 = sqrt(dx2^2 + dy2^2),
                     n = n()
  )) %>% 
  group_by(id) %>% 
  mutate(index = seq(1, n(), 1))

## KS test
collapsed <- filter(lag4, !is.na(dt4)) %>%
  mutate(ToD2 = forcats::fct_collapse(ToD, dawn = "dawn", group_other = T))
classes <- collapsed %>% dplyr::select(ToD2, dt4) %$% downSample(dt4, ToD2)
dawn_class <- filter(classes, Class == "dawn") %>% pull(x)
other_class <- filter(classes, Class == "Other") %>% pull(x)
ks.test(dawn_class, other_class)

## Time of day points near water?
buffers <- st_buffer(water_utm, 500)
water_pts <- st_intersection(buffers, 
                             st_as_sf(lag4, 
                                      coords = c("x", "y"), 
                                      crs = 32733, na.fail = F))
water_pts %>% select(-n) %>% group_by(ToD) %>% tally()

