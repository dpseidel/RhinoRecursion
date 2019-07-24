### Code for 24 hour displacement analysis associated with Seidel et al. 2019 
### "Mesoscale movement and recursion behaviors of Namibian black rhinos"
### (submitted to Movement Ecology, BioMove Special Edition)
### Author: Dana Paige Seidel

source("DataParsing_Script.R")
library(caret)

## Time of Day displacement analysis
sampling_days <- c(seq.Date(from = ymd("2011-10-19"), to = ymd("2014-01-24"), 1), 
                   seq.Date(from = ymd("2017-06-20"), to = ymd("2018-11-29"), 1))
intervals <- interval(start = ymd_hm(paste(sampling_days - 1, "22:30")), 
                      end = ymd_hm(paste(sampling_days, "22:30")))

ToDfilter <- all_rhinos %>% st_set_geometry(NULL) %>% 
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
         offset_hour = case_when(
           rnd_hour == 0 & hour(date) == 23 ~ difftime(ymd_hms(paste(local_date, local_time)), 
                                                       ymd_h(paste(local_date + 1, rnd_hour)), 
                                                       units = "secs"), 
           TRUE ~ difftime(ymd_hms(paste(local_date, local_time)),
                           ymd_h(paste(local_date, rnd_hour)),
                           units = "secs") 
         ),
         offset_ToD = case_when(
           rnd_hour %in% c(0, 23) & hour(date) != 0 ~ difftime(ymd_hms(paste(local_date, local_time)), 
                                                               ymd_h(paste(local_date + 1, ToD_est)), 
                                                               units = "secs"), 
           TRUE ~ difftime(ymd_hms(paste(local_date, local_time)),
                           ymd_h(paste(local_date, ToD_est)),
                           units = "secs") 
         )) %>% 
  filter(ToD != "other") %>% 
  mutate(interval = map_dbl(date, ~which(.x %within% intervals))) %>% 
  group_by(id, interval, ToD) %>%  
  mutate(closest =  abs(as.numeric(offset_ToD)) == min(abs(as.numeric(offset_ToD)))) %>% 
  filter(closest) %>% ungroup() %>% 
  mutate(ToD = readr::parse_factor(ToD, ordered = T, levels = c("mid-night", "dawn", "mid-day", "dusk")))

ToD <- ToDfilter %>% split(.$id) %>% 
  map(., ~.x %>%
        complete(interval = seq(from = min(interval, na.rm = T), to = max(interval, na.rm = T), 1)) %>% 
        complete(ToD, nesting(interval)) %>% filter(!is.na(ToD)) %>% 
        fill(id) %>% 
        arrange(interval, ToD)) 

# test that the ordering happened correctly, crucial to next step:
sum(map_lgl(ToD,~ all(.x$ToD == rep_len(levels(ToDfilter$ToD), length(.x$ToD))))) == 59

lag4 <- ToD %>%
  map_df(., function(x) {
    
    df <- mutate(x,
           dx = c(diff(x), NA),
           dy = c(diff(y), NA),
           dt = sqrt(dx^2 + dy^2),
           dx2 = c(diff(x, lag = 2), NA, NA),
           dy2 = c(diff(y, lag = 2), NA, NA),
           dt2 = sqrt(dx2^2 + dy2^2),
           dx4 = c(diff(x, lag = 4), NA, NA, NA, NA),
           dy4 = c(diff(y, lag = 4), NA, NA, NA, NA),
           dt4 = sqrt(dx4^2 + dy4^2)
    )
    
    # extracting true time between points
    # being explicit about units is important. 
    s <- diff(df$date)
    s2 <- diff(df$date, lag = 2)
    s4 <- diff(df$date, lag = 4)  
    units(s) <- "hours"
    units(s2) <- "hours"
    units(s4) <- "hours"
    
    df$time <- c(as.vector(s), NA)
    df$time2 <- c(as.vector(s2), NA, NA)
    df$time4 <- c(as.vector(s4), NA, NA, NA, NA)
    
    return(df)
  })

### make sure time diffferences are as expected:
mean(lag4$time, na.rm = T)
range(lag4$time, na.rm = T)
sd(lag4$time, na.rm = T)

mean(lag4$time2, na.rm = T)
sd(lag4$time2, na.rm = T)
range(lag4$time2, na.rm = T)

mean(lag4$time4, na.rm = T)
range(lag4$time4, na.rm = T)
sd(lag4$time4, na.rm = T)

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

### Periodicity test!
chosen9 %>% split(.$id) %>% 
  map_dbl(nrow)
# a couple have 105-155 relocations but most < 100

chosen9 %>% split(.$id) %>% 
  map(function(x){ptest::ptestg(na.omit(x$dt2))})
# rejected the null hypothesis of periodicity (@ alpha = .05) for 1 individual SAT2372.  

## KS test
collapsed <- filter(lag4, !is.na(dt4)) %>%
  mutate(ToD2 = forcats::fct_collapse(ToD, dawn = "dawn", group_other = T))
classes <- collapsed %>% dplyr::select(ToD2, dt4) %$% downSample(dt4, ToD2)
dawn_class <- filter(classes, Class == "dawn") %>% pull(x)
other_class <- filter(classes, Class == "Other") %>% pull(x)
ks.test(dawn_class, other_class)

## Time of day points near water?
buffers <- st_buffer(water_utm, 250)
water_pts <- st_intersection(buffers, 
                             st_as_sf(lag4, 
                                      coords = c("x", "y"), 
                                      crs = 32733, na.fail = F))
water_pts %>% group_by(ToD) %>% tally()

