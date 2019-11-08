Cairo::CairoPNG(filename = "stddisp.png", width = 1200, height = 1000,
                res = 180)
ggplot(lag4) +
  geom_density(aes(x = dt4/time4, color = ToD, fill = ToD), alpha = .1, na.rm = T) +
  theme_minimal() +
  theme(legend.position = "top", plot.title = element_text(hjust = .5)) +
  labs(title = "Displacement per hour distributions across time of day",
       x = "24-hr displacement per hour (m/hr)") +
  scale_color_manual("Time of Day", values = c("#26384f", "#57beff", "#4f8a99", "#00578b")) +
  scale_fill_manual("Time of Day", values = c("#26384f", "#57beff", "#4f8a99", "#00578b"))
dev.off()

