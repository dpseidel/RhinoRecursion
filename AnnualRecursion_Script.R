### Code for Annual recursion analysis associated with Seidel et al. 2019 
### "Mesoscale movement and recursion behaviors of Namibian black rhinos"
### (submitted to Movement Ecology, BioMove Special Edition)
### Author: Dana Paige Seidel

source("DataParsing_Script.R")
library(tlocoh)
library(tlocoh.dev)

year_rhinos <- all_rhinos %>% group_by(id) %>% summarize(
  start = min(date), end = max(date), 
  days = difftime(end, start)
) %>% filter(days > 360) %>% pull(id)

rec_rhinos <- all_rhinos %>% 
  filter(id %in% year_rhinos[c(1:2, 4:7)],
         date >= ymd("2012-04-01"), 
         date < ymd("2013-04-01")) %>% 
  st_set_geometry(NULL) %>% split(.$id)

map_tu <- function(df, ivg){
  df <- na.omit(df)
  k <-  round(sqrt(nrow(df)))
  lxy <- xyt.lxy(xy = matrix(c(df$x, df$y), ncol = 2), 
                 dt = df$date,
                 id = df$id, 
                 proj4string = CRS("+init=epsg:32733"))
  lxy.tumap(lxy, ivg=ivg, grid = "square", cellsize = 1000)
}

tumaps_7day <- map(rec_rhinos, map_tu, ivg = 7*24*3600)
tumaps_sf <- map(tumaps_7day, ~.x[[1]] %>% st_as_sf)

# Look at full set of time use grids, identify/drop watering holes
grids <- pmap(list(x = tumaps_sf, y= names(tumaps_sf)), 
              function(x, y){
                t <- x %>% 
                  mutate(waterhole = as.numeric(st_intersects(., water_utm, sparse = T)), 
                         grid.id = 1:nrow(.), id = y) %>% 
                  filter(nsv.604800 > 2, is.na(waterhole))
                return(t)
              })

rhino_list <- all_rhinos %>% filter(id %in% year_rhinos[1:7]) %>% split(.$id)

joins <- pmap(list(grid = grids, pts = rhino_list[c(1:2, 4:7)]), 
              function(grid, pts){
                st_join(pts, st_transform(grid, st_crs(pts))) %>% filter(!is.na(grid.id))
              })

df_3visits <- map_df(joins,~ group_by(.x, grid.id) %>% mutate(
  diff = c(NA, as.numeric(diff(date))),
  units = units(diff(date)),
  diff_days = case_when(
    units == "days" ~ diff, 
    units == "hours" ~ diff/24 
  ), 
  sv = ifelse(is.na(diff), T, diff_days > 7),
  nsv = sum(sv)) %>% 
    st_set_geometry(NULL))


### Correlation Tests
# grid sizes. 
map_dbl(tumaps_sf, ~nrow(.x))

# but some of those extents include a huge amount of area not used, not even walked in. 
visited <- map_dbl(tumaps_sf, ~filter(.x, nsv.604800 > 0) %>% nrow())

# revisited
revisited <- map_dbl(tumaps_sf, ~filter(.x, nsv.604800 >= 2) %>% nrow())

# proportion revisited
revis_prop <- revisited/visited

# mean number of returns to revisited cells. 
avg_revis <- map_dbl(tumaps_sf, ~filter(.x, nsv.604800 >= 2) %$% mean(nsv.604800, na.rm = T))
sd_revis <- map_dbl(tumaps_sf, ~filter(.x, nsv.604800 >= 2) %$% sd(nsv.604800, na.rm = T))
med_revis <- map_dbl(tumaps_sf, ~filter(.x, nsv.604800 >= 2) %$% median(nsv.604800, na.rm = T)) 

tumap_stats <- tibble(id = names(revis_prop), revis_prop, visited, avg_revis, sd_revis, med_revis)

cor(tumap_stats[,-1])

cor(revis_prop, visited)
cor(revis_prop[-6], visited[-6])  # 280 is an outlier... not sure why. 

cor(med_revis, visited)

# range size has a negative relationship with proportion revisted cells
# i.e. smaller ranges are likely to revisited more of their range. 
# our analysis is of course, limited by sample size. though this is a logical result. 

# number of cells within the top 25% of revisits. 
top25rev <- map_dbl(tumaps_sf, ~filter(.x, nsv.604800 > 0)  %>% 
                      filter(., nsv.604800 > quantile(nsv.604800, probs = .75)) %>% 
                      nrow())

# the proportion of cells visited that gain the top 25% of returns, 
# is fairly equal across individuals, and is not strongly correlated with range size. 
cor(top25rev, visited)