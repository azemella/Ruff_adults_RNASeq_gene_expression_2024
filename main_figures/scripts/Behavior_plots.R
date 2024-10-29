# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# Read data files for aggression and courtship morphs
agg <- read.table("data/hormones_and_behaviors_data/Aggression_morphs_RNAseq.txt", header = TRUE, sep = "\t")
court <- read.table("data/hormones_and_behaviors_data/Courtship_morphs_RNAseq.txt", header = TRUE, sep = "\t")

# Convert the 'morph' column in both datasets to factors
agg$morph <- as.factor(agg$morph)
court$morph <- as.factor(court$morph)

# Reorder factor levels
agg$morph <- factor(agg$morph, levels = c("Ind", "Sat", "Fae"))
court$morph <- factor(court$morph, levels = c("Ind", "Sat", "Fae"))

# Create ggplot object
plotagg <- agg %>%
  ggplot(aes(y = aggressionrate, x = morph)) + 
  geom_point(size = -1, aes(fill = morph, color = morph)) +
  geom_boxplot(aes(fill = morph), 
               position = position_dodge(0.9), 
               width = 0.8, 
               size = 0.35,
               linewidth = 0.5,
               outlier.shape = NA,
               show.legend = FALSE) +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper.., fill = morph), 
               size = 0.35,
               width = 0.8,
               linewidth = 0.5,
               position = position_dodge(0.9),
               show.legend = FALSE,
               outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax.., fill = morph), linewidth = 0.5, width = 0.5, size = 0.35, alpha = 1, 
               position = position_dodge(0.9), show.legend = FALSE, outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin.., fill = morph), linewidth = 0.5, width = 0.5, size = 0.35, alpha = 1, 
               position = position_dodge(0.9), show.legend = FALSE, outlier.shape = NA) +
  geom_point(aes(fill = morph), 
             size = 2, 
             shape = 21, 
             position = position_jitterdodge(jitter.width = 0.5, 
                                             dodge.width = 0.9)) +
  geom_signif(data = data,
              aes(y_position = 8.25, 
                  xmin = 0.75, 
                  xmax = 2.25,
                  annotations = "***"), 
              tip_length = 0.01, 
              manual = T,
              vjust = 0.65,
              size = 0.35,
              textsize = 7,
              inherit.aes = TRUE) +
  geom_signif(data = data,
              aes(y_position = 8.75, 
                  xmin = 0.75, 
                  xmax = 3.2,
                  annotations = "***"), 
              tip_length = 0.01, 
              manual = T,
              vjust = 0.65,
              size = 0.35,
              textsize = 7,
              inherit.aes = TRUE) +
  scale_fill_manual(labels = c("Independent","Satellite","Faeder"),
                    values = c("blue", "violetred", "orange")) +
  scale_color_manual(labels = c("Independent","Satellite","Faeder"),
                     values = c("blue", "violetred", "orange")) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8), 
                     limits = c(0, 8.75), 
                     labels = c("0", "2", "4", "6", "8")) +
  guides(color = guide_legend(override.aes = list(shape = 21, size = 6))) +
  labs(x = "",
       y = expression("Aggressive behaviors min"^{-1})) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 11.5, colour = "black", family = "Arial"),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 18, family = "Arial", margin = margin(t = 0, r = 8, b = 0, l = 0)),
    axis.title.x = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.75),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = alpha("white", 1)),
    legend.title = element_blank(),
    legend.text = element_text(size = 18, family = "Arial"),
    legend.position = "bottom",
    legend.key = element_blank(),
    legend.box.just = "left",
    legend.background = element_blank(),
    legend.direction = "horizontal",
    panel.spacing = unit(0.3, "lines"),
    strip.background = element_rect(fill = "white"),
    strip.text.x = element_text(size = 20, family = "Arial"))

# Create ggplot object
plotcourt <- court %>%
  ggplot(aes(y = courtshiprate, x = morph)) + 
  geom_point(size = -1, aes(fill = morph, color = morph)) +
  geom_boxplot(aes(fill = morph), 
               position = position_dodge(0.9), 
               width = 0.8, 
               size = 0.35,
               linewidth = 0.5,
               outlier.shape = NA,
               show.legend = FALSE) +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper.., fill = morph), 
               size = 0.35,
               width = 0.8,
               linewidth = 0.5,
               position = position_dodge(0.9),
               show.legend = FALSE,
               outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax.., fill = morph), linewidth = 0.5, width = 0.5, size = 0.35, alpha = 1, 
               position = position_dodge(0.9), show.legend = FALSE, outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin.., fill = morph), linewidth = 0.5, width = 0.5, size = 0.35, alpha = 1, 
               position = position_dodge(0.9), show.legend = FALSE, outlier.shape = NA) +
  geom_point(aes(fill = morph), 
             size = 2, 
             shape = 21, 
             position = position_jitterdodge(jitter.width = 0.5, 
                                             dodge.width = 0.9)) +
  geom_signif(data = data,
              aes(y_position = 10, 
                  xmin = 0.85, 
                  xmax = 2.25,
                  annotations = "*"), 
              tip_length = 0.01, 
              manual = T,
              vjust = 0.65,
              size = 0.35,
              textsize = 7,
              inherit.aes = TRUE) +
  geom_signif(data = data,
              aes(y_position = 10.5,
                  xmin = 0.85, 
                  xmax = 3.2,
                  annotations = "**"), 
              tip_length = 0.01, 
              manual = T,
              vjust = 0.65,
              size = 0.35,
              textsize = 7,
              inherit.aes = TRUE) +
  geom_signif(data = data,
              aes(y_position = 9.5, 
                  xmin = 2.25, 
                  xmax = 3.2,
                  annotations = "**"), 
              tip_length = 0.01, 
              manual = T,
              vjust = 0.65,
              size = 0.35,
              textsize = 7,
              inherit.aes = TRUE) +
  scale_fill_manual(labels = c("Independent","Satellite","Faeder"),
                    values = c("blue", "violetred", "orange")) +
  scale_color_manual(labels = c("Independent","Satellite","Faeder"),
                     values = c("blue", "violetred", "orange")) +
  guides(color = guide_legend(override.aes = list(shape = 21, size = 6))) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10), 
                     limits = c(0, 10.5), 
                     labels = c("0", "2", "4", "6", "8", "10")) +
  labs(x = "",
       y = expression("Courtship behaviors min"^{-1})) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 11.5, colour = "black", family = "Arial"),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 18, family = "Arial", margin = margin(t = 0, r = 8, b = 0, l = 0)),
    axis.title.x = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.75),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = alpha("white", 1)),
    legend.title = element_blank(),
    legend.text = element_text(size = 18, family = "Arial"),
    legend.position = "bottom",
    legend.key = element_blank(),
    legend.box.just = "left",
    legend.direction = "horizontal",
    legend.background = element_blank(),
    panel.spacing = unit(0.3, "lines"),
    strip.background = element_rect(fill = "white"),
    strip.text.x = element_text(size = 20, family = "Arial"))

# Save plot as SVG
arrange <- ggarrange(plotTT, plotA4, plotagg, plotcourt, nrow = 1, ncol = 4, common.legend = TRUE, legend = "bottom")
ggsave(filename = "main_figures/plots/Hormones_and_behaviours.svg", plot = arrange, width = 8, height = 4, 
       bg = "white", dpi = "print", device = "svg")
