### Code for figures associated with Seidel et al. 2019
### "Mesoscale movement and recursion behaviors of Namibian black rhinos"
### (submitted to Movement Ecology, BioMove Special Edition)
### Author: Dana Paige Seidel

####### Necessary Libraries #######
library(stmove) # package not available on CRAN, see github.com/dpseidel/stmove
library(ggsn) # for spatial scale and coordinate arrow

########### Source Analysis Scripts ##########
source("24hrDisplacement_Script.R") # defines: chosen9, lag4, ToD
source("16dayRecursion_HomeRange_Script.R") # defines: biweekly, tumaps
source("AnnualRecursion_Script.R") #  defines: df_3visits

####### Figures #######

#### Figure 1
means_rhinos <- stmove::dist_map(
  st_set_geometry(all_rhinos, NULL),
  "+proj=utm +zone=33 +south"
)
ggplot() +
  geom_sf(data = namib_eco, aes(fill = ECO_NAME)) +
  scale_fill_manual("", # get everyone rearranged in a way that makes sense
    values = terrain.colors(8)[c(1, 8, 2, 4, 5:7, 3)],
    guide = guide_legend(override.aes = list(shape = NA))
  ) +
  ggsn::scalebar(
    data = namib_eco, dist = 100, dist_unit = "km",
    transform = TRUE, model = "WGS84",
    anchor = c(x = 16.5, y = -21.25), st.dist = .05
  ) +
  ggsn::north(data = namib_eco, anchor = c(x = 13.5, y = -20.6), symbol = 16, scale = .2) +
  geom_sf(data = means_rhinos, size = 1, aes(shape = "shape"), show.legend = "point") +
  scale_shape_manual("", labels = "Mean Rhino Location", values = 19) +
  theme_void() #+ 
# labs(caption = "Ecoregions as defined by Dinerstein et al. 2017") +
# theme(plot.caption = element_text(hjust = 2))

#### Figure 2 - Time Line
stmove::plot_timeline(all_rhinos)

#### Figure 3 - Dawn displacement time series
Cairo::CairoPNG(filename = "Figures/distplacement.png", width = 1500, height = 1500,
                res = 180)
chosen9 %>%
  ungroup() %>%
  filter(ToD == "dawn") %>%
  mutate(
    id_f = forcats::fct_reorder(id, yday(local_date)),
    sex = case_when(
      id %in% paste0("SAT", c("2369", "2372", "2414", "278", "680")) ~ "F",
      id %in% paste0("SAT", c("277", "280", "682", "2058")) ~ "M"
    )
  ) %>%
  group_by(id_f) %>%
  mutate(
    x_lab = min(yday(date)),
    nudge = sd(yday(date))
  ) %>%
  ggplot() +
  geom_path(aes(disphr, x = yday(date), color = sex), na.rm = T) +
  scale_color_manual(values = c("black", "#000080")) +
  # geom_text(aes(x = x_lab + nudge/3, y = 15500,
  #              label = paste(round(n/2), "days")), size = 3) +
  coord_cartesian(ylim = c(0, 700)) +
  facet_wrap(~id_f, scales = "free_x", ncol = 3) +
  theme_minimal() +
  labs(
    title = "Dawn-Dawn standardized 24-hr displacement through time",
    x = "Julian day", y = "24-hr displacement per hour (m/hr)",
    caption = "Note: A single observation approaching 1250 m/hr displacement was cut off from SAT277's graph to maintain comparable scales"
  ) +
  theme(plot.title = element_text(hjust = .5), legend.position = "none", panel.spacing = unit(1, "lines"))
dev.off()

#### Figure 4. Facet Distributions -- distributions
Cairo::CairoPNG(filename = "Figures/distributions.png", width = 1500, height = 1500,
                res = 180)

chosen9 %>%
  ungroup() %>%
  filter(ToD == "dawn") %>%
  mutate(
    id_f = forcats::fct_reorder(id, yday(date)),
    sex = case_when(
      id %in% paste0("SAT", c("2369", "2372", "2414", "278", "680", "676", "677")) ~ "F",
      id %in% paste0("SAT", c("277", "280", "678", "682", "2058")) ~ "M"
    )
  ) %>%
  group_by(id_f) %>%
  mutate(
    x_lab = min(yday(date)),
    nudge = sd(yday(date))
  ) %>%
  ggplot() +
  geom_density(aes(disphr, color = sex), na.rm = T) +
  scale_color_manual(values = c("black", "#000080")) +
  geom_text(aes(
    x = 900, y = .0125,
    label = paste(round(n / 2), "days")
  ), size = 3) +
  facet_wrap(~id_f, ncol = 3) +
  theme_minimal() +
  labs(
    title = "Dawn-Dawn standardized 24-hr displacement through time",
    x = "24-hr displacement per hour (m/hr)"
  ) +
  theme(plot.title = element_text(hjust = .5), legend.position = "none", panel.spacing = unit(1, "lines"),
        plot.margin = margin(10, 15, 10, 10))
dev.off()


####  Fig. 5
Cairo::CairoPNG(filename = "Figures/density.png", width = 1500, height = 1500,
                res = 180)
ggplot(lag4) +
  geom_density(aes(x = dt4/time4, color = ToD, fill = ToD), alpha = .1, na.rm = T) +
  theme_minimal() +
  theme(legend.position = "top", plot.title = element_text(hjust = .5)) +
  labs(title = "Standardized 24-hr displacement distributions across time of day", 
       x = "24-hr displacement per hour (m/hr)") +
  scale_color_manual("Time of Day", values = c("#26384f", "#57beff", "#4f8a99", "#00578b")) +
  scale_fill_manual("Time of Day", values = c("#26384f", "#57beff", "#4f8a99", "#00578b"))
dev.off()

#### Figure 6.
# T-Locoh  TimeUse Plot.
par(mfrow = c(1, 2))
plot(tumaps$`SAT189-2011-305`$result)

#### Figure 7
filter(biweekly, nsv.43200 > 1) %>%
  ggplot(aes(y = nsv.43200, x = meanNDVI * .0001)) +
  stat_bin2d(aes(fill = log(stat(count))), binwidth = c(.05, 1)) +
  # stat_bin2d(aes(color = ..count..), size = 1, bins = 30) +
  scale_fill_gradient("Count",
    low = "#c6dbef", high = "#08306b",
    labels = function(x) {
      round(exp(1)^x)
    }
  ) +
  # geom_jitter(size = .6, alpha = .3) +
  # scale_color_viridis_c(option="inferno")+
  labs(
    y = "Number of Separate Visits",
    x = "Average NDVI",
    title = "Relationship between 16-day revisitation and grid cell NDVI"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = .5))
dev.off()

#### Figure 8:
ggplot(filter(df_3visits, sv == T) %>% group_by(id.x) %>% 
         mutate(num = length(unique(grid.id)))) + 
  geom_histogram(aes(x = diff_days, fill = cut(nsv, breaks = 10, dig.lab = 1)), binwidth = 7) + 
  geom_text(x= 275, y = 105, aes(label = paste(num, "cells")), size = 3.5) +
  facet_wrap(. ~ as.factor(id.x)) +
  scale_fill_manual("Number of\nseparate visits", 
                    values = colorRampPalette(c("#08306b", "#9ecae1"))(9)) +
  theme_minimal() +
  labs(x = "Time between visits (days)", 
       title = "Annual recursion: time to return across individuals") +
  theme(plot.title = element_text(hjust = .5))
