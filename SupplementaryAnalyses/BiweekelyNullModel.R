##### Null Distribution of NDVI -- for 16 day revisitation grids?

### randomly sample -- the same number of points from each NDVI file
### and compare to the biweekly extracted points. 

### do we want to randomly sample within a home range?
set.seed(4236)

biweekly_null <- biweekly %>% 
  group_by(id, year, date_id) %>% 
  filter(nsv.43200 > 0, !is.na(meanNDVI)) %>% 
  mutate(nullNDVI_idint = sample(meanNDVI, replace = F)) %>% 
  group_by(id) %>% 
  mutate(nullNDVI_id = sample(meanNDVI, replace = F)) %>% 
  ungroup() %>% 
  mutate(nullNDVI = sample(meanNDVI, replace = F))


Cairo::CairoPNG(filename = "ndvidist.png", width = 1200, height = 1000,
                res = 180)
ggplot() + 
  geom_histogram(aes(x = meanNDVI * .0001, 
                     fill = ifelse(nsv.43200 > 0, "visited", "not-visited" )), 
                 data = biweekly, bins = 100) + 
  scale_fill_manual("Visitation", values = c("#696969", "#4169e1")) +
  theme_minimal() + labs(
    title = "Distribution of extracted NDVI values", 
    x = "Average NDVI within extracted 1km^2 grid cell"
  )  +
  theme(plot.title = element_text(hjust = .5))
dev.off()

library(patchwork)
p_org = filter(biweekly_null, nsv.43200 > 1) %>%
  ggplot(aes(y = nsv.43200, x = meanNDVI * .0001)) +
   stat_bin2d(aes(fill = log(stat(count))), binwidth = c(.05, 1)) +
  # # stat_bin2d(aes(color = ..count..), size = 1, bins = 30) +
   scale_fill_gradient("Count",
                       low = "#c6dbef", high = "#08306b",
                       labels = function(x) {
                         round(exp(1)^x)
                       }
   ) +
  # geom_jitter(size = .6, alpha = .3) +
  #geom_smooth(method = "glm") +
  #scale_x_continuous(limits = c(-0.149, 0.725)) +
  # scale_color_viridis_c(option="inferno")+
  labs(
    y = "Number of Separate Visits",
    x = "Average NDVI",
    title = "Null Model -- all visited cells"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = .5))

filter(biweekly_null, nsv.43200 > 1) %>%
  ggplot(aes(y = nsv.43200, x = nullNDVI_id * .0001)) +
  stat_bin2d(aes(fill = log(stat(count))), binwidth = c(.05, 1)) +
  # stat_bin2d(aes(color = ..count..), size = 1, bins = 30) +
  scale_fill_gradient("Count",
                      low = "#c6dbef", high = "#08306b",
                      labels = function(x) {
                        round(exp(1)^x)
                      }
  ) +
  # geom_jitter(size = .6, alpha = .3) +
  #geom_smooth(method = "glm", color = "red") +
  # scale_color_viridis_c(option="inferno")+
  labs(
    y = "Number of Separate Visits",
    x = "Average NDVI",
    title = "Null Model -- ID specific"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = .5))


filter(biweekly_null, nsv.43200 > 1) %>%
  ggplot(aes(y = nsv.43200, x = nullNDVI_id * .0001)) +
   geom_jitter(size = .6, alpha = .3) +
  stat_smooth(method = "glm", formula = y ~ poly(x, 2), size = 1) +
  geom_smooth(method = "glm") +
  # scale_color_viridis_c(option="inferno")+
  labs(
    y = "Number of Separate Visits",
    x = "Average NDVI",
    title = "Null Model - per individual"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = .5))

lm <- lm(`nsv.43200` ~ (meanNDVI*.0001),  data = biweekly_null)
quad <- lm(nsv.43200 ~ poly(meanNDVI*.0001, 2),  data = biweekly_null)

anova(lm, quad)

filter(biweekly_nullrev, nsv.43200 > 1) %>%
  ggplot(aes(y = nsv.43200, x = meanNDVI * .0001)) +
  geom_jitter(size = .6, alpha = .3) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 2), size = 1)
  #geom_smooth(method = "glm") +
  #geom_smooth(aes(y = nsv.43200, x = nullNDVI * .0001), method = "glm", col = "red") +
  # scale_color_viridis_c(option="inferno")+
  labs(
    y = "Number of Separate Visits",
    x = "Average NDVI",
    title = "Null Model - per individual per interval"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = .5))



biweekly_nullrev %>% 
  filter(nsv.43200 > 1) %>% 
  summarise(avg = mean(meanNDVI, na.rm = T), sd = sd(meanNDVI, na.rm = T),
            null = mean(nullNDVI, na.rm = T), nsd = sd(nullNDVI, na.rm = T)) %>% 
clipr::write_clip()

biweekly_nullid %>% ungroup() %>%  
  filter(nsv.43200 > 1) %>% 
  summarise(avg = mean(meanNDVI, na.rm = T), sd = sd(meanNDVI, na.rm = T),
            null = mean(nullNDVI, na.rm = T), nsd = sd(nullNDVI, na.rm = T)) %>% 
  clipr::write_clip()


### If you want to find an intermediate value you have to test the fit of a
### quadratic y=ax^2 + bx + c and see if this fit better from the AIC point of 
### view than your linear fit.  

mod_biweekly <- biweekly_null %>% 
  mutate(meanNDVI = .0001*meanNDVI, 
         nullNDVI_id = .0001*nullNDVI_id
         )

mod1 <- glm(nsv.43200 ~ nullNDVI_id, data = mod_biweekly, family = "poisson") 
mod2 <- glm(nsv.43200 ~ meanNDVI, data = mod_biweekly, family = "poisson") 
mod3 <- glm(nsv.43200 ~ poly(meanNDVI, 2), data = mod_biweekly, family = "poisson") 


anova(mod1, mod2, mod3)
AIC(mod1, mod2, mod3)

Cairo::CairoPNG(filename = "quad.png", width = 1200, height = 1000,
                res = 180)
mod_biweekly %>%
  ggplot(aes(y = nsv.43200, x = meanNDVI * .0001)) +
  geom_jitter(size = .6, alpha = .3) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 2), size = 1) +
  labs(
  y = "Number of Separate Visits",
  x = "Average NDVI",
  title = "Average NDVI values in visited cells, fitted with a Quadratic Model"
) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = .5))
dev.off()
