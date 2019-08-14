##### Group effects, Seasonality


### Seasonal Differences ###
# wet season defined as november - april, according to rainfall data
library(readxl)
rainfall1981_2013 <- read_xlsx("Data/ENP_rainfall_daily_1981-2013_clean.xlsx",
                               col_types = c(
                                 "date", "numeric", "text", "numeric", # ignore the parsing errors
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"
                               )
)

Cairo::CairoPNG(filename = "rainfall.png", width = 1200, height = 1000,
                res = 180)
rainfall1981_2013 %>%
  gather(., area, rainfall, -month, -year, -month_cal, -day, -date) %>%
  dplyr::group_by(month, year, area) %>%
  dplyr::summarise(rainfall_month = sum(rainfall, na.rm = T)) %>%
  group_by(month) %>%
  summarise(mean_rainfall = mean(rainfall_month, na.rm = T)) %>%
  ggplot(., aes(x = month, y = mean_rainfall)) + geom_bar(stat = "identity") +
  ggtitle("Average Monthly Rainfall across all areas ENP, 1981-2013") +
  theme_minimal() +
  scale_x_continuous(breaks = c(1:12), 
                     labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  labs(y = "Average Monthly Rainfall (mm)") +
  theme(plot.title = element_text(hjust = .5))
dev.off()

# wet == November - April 
wet_chosen <- chosen9 %>% mutate(wet = ifelse(month(local_date) %in% c(1:4, 11:12), T, F))
wet_lag <- lag4 %>% mutate(wet = ifelse(month(local_date) %in% c(1:4, 11:12), T, F))

Cairo::CairoPNG(filename = "season9.png", width = 1200, height = 1000,
                res = 180)
ggplot(wet_chosen, aes(x = dt2, fill = wet)) +
  geom_density(na.rm = T, alpha = .4) + 
  labs(title = "24-hr Displacement Across Seasons - 9 Consistent Individuals", 
       x = "24-hr displacement (m)") + 
  scale_fill_discrete("Season", labels = c("dry", "wet")) +
  facet_wrap(~ToD) +
  theme_minimal() + 
  theme(strip.background = element_rect(fill = "#F4F4F4"),
        plot.title = element_text(hjust = .5))
dev.off()

wet_chosen %>% select(-n) %>% group_by(wet) %>% tally()

Cairo::CairoPNG(filename = "season_TOD.png", width = 1200, height = 1000,
                res = 180)
ggplot(wet_lag, aes(x = dt4, fill = wet)) + 
  geom_density(alpha = .4, na.rm = T) + facet_wrap(~ToD) +
  theme_minimal() + 
  labs(title = "24-hr Displacement Across Seasons", 
       x = "24-hr displacement (m)") +
  scale_fill_discrete("Season", labels = c("dry", "wet")) +
  theme(strip.background = element_rect(fill = "#F4F4F4"),
        plot.title = element_text(hjust = .5))
dev.off()

wet_lag %>% na.omit() %>% group_by(wet) %>% tally() 

#### Demographics ####
ids <- ToDfilter %>% pull(id) %>% unique()
metadataDbbCss <- readxl::read_xlsx("../../../Dropbox/Rhino/Collar Info/Black rhino sat collars/Rhino Sat Collar Data/RhinoSAT_combined_editted.xlsx", sheet = 1)

metadataDbb1718 <- readxl::read_xlsx("../../../Dropbox/Rhino/Rhino Movements/moredataetosha/ENPDbbmove2017_18.xlsx", sheet = 1)

# bring it all together, sort out calf presence/pregnant. 
old_collars <- metadataDbbCss %>% 
  mutate(id = paste0(Type, `Collar ID`), 
         id = ifelse(id == "SAT843", "SAT943", id)) %>%  # I have a hunch collar "943" is actually listed as collar "843" in this metadata. patching.. 
  filter(id %in% ids) 

Dbb1718_collars <- metadataDbb1718 %>% 
  mutate(id = str_remove(str_remove(`Tag Number`, "IR-"), " ")) %>% 
  filter(id %in% ids)

nrow(old_collars) + nrow(Dbb1718_collars) == length(ids)
ids %in% c(old_collars$id, Dbb1718_collars$id)

collars_demographics <- old_collars %>% 
  select(id, Sex, Composition, Age, Remarks) %>%
  bind_rows(., Dbb1718_collars %>% 
              select(id, Sex = Gender, Remarks = Notes)) %>% 
  mutate(has_calf = case_when(
    str_detect(Composition, "Size") ~ TRUE, 
    str_detect(Composition, "calf") ~ TRUE, 
    str_detect(Remarks,  "calf") ~ TRUE, 
    TRUE ~ FALSE
  ), 
  pregnant = str_detect(Remarks, "pregnant"), 
  pregnant = ifelse(is.na(pregnant), FALSE, pregnant), 
  Sex = case_when(
    str_detect(Sex, "^F") ~ "F", 
    str_detect(Sex, "^M") ~ "M", 
    TRUE ~ NA_character_
  ), 
  Condition = case_when(
    has_calf ~  "with calf", 
    pregnant ~ "with calf", 
    !has_calf ~ "single")) 


# When pregant females and females captured with calf are combine we have nearly even
# groups. Note that one individual's, SAT2596, sex was not recorded and was dropped from
# the analysis. 

collars_demographics %>% filter(!is.na(Sex)) %>% 
  group_by(Sex, Condition) %>% 
  tally() %>% arrange(Sex) %>% knitr::kable("latex")
  
demo_lag <- lag4 %>% right_join(collars_demographics) %>% filter(!is.na(Sex))

Cairo::CairoPNG(filename = "sexTOD.png", width = 1200, height = 1000,
                res = 180)
ggplot(demo_lag, aes(x = dt4, fill = ToD)) +
  geom_density(na.rm = T, alpha = .4) + facet_wrap(~Condition +Sex) +
  theme_minimal() + 
  labs(x = "24-hour displacement") +
  theme(strip.background = element_rect(fill = "#F4F4F4"),
        plot.title = element_text(hjust = .5))
dev.off()

ggplot(filter(demo_lag, str_detect(ToD, "^d")), aes(x = dt4, fill = ToD)) +
  geom_density(na.rm = T, alpha = .4) + facet_wrap(~Condition +Sex)



##### Regional Differences

# Kunene

# Etosha --- across west/east/central 

# Waterberg --  only 4 individuals. 
library(stmove)
means <- dist_map(st_set_geometry(all_rhinos, NULL), proj = "+proj=utm +zone=33 +south +datum=WGS84 +units=m +no_defs")
park <- read_sf("../../../Dropbox/Research - PHD/Anthrax/ENP data for transmission study/ENP shapefiles/ENP_three_areas.shp") %>% 
  st_transform(32733) %>% mutate(id = c("East", "West", "Central"))

groups <- st_intersects(means, park, sparse = F) %>% 
  as_tibble(.name_repair = "unique") %>% 
  rename(East = "...1", West = "...2", Central = "...3") %>%
  mutate(ids = means$id, 
         Area = case_when(
           East == TRUE ~ "East",
           West == TRUE ~ "West", 
           Central == TRUE  ~ "Central", 
           TRUE ~ "OOP"
         ))

waterberg_ids <- paste0("SAT", c(469, 676:678, 680:683))
groups[groups$ids %in% waterberg_ids, "Area"] <- "Waterberg"

groups %>% group_by(Area) %>% tally()

Rhino_meanlocs <- left_join(means, groups, by = c("id" = "ids"))
