# Figure 2: saved-model predictions.

library(lme4)
library(ggeffects)
library(emmeans)
library(ggplot2)
library(ggsci)
library(patchwork)

# Set the working directory to the repository root before running.
dir.create("results", showWarnings = FALSE)
best_model <- readRDS("models/Step4_best_model_Results.rds")
df_agg <- model.frame(best_model)

model_weights <- rowSums(model.response(df_agg))
ref_nat <- weighted.mean(df_agg$C1_NatStatus, model_weights)
ref_invasive <- weighted.mean(df_agg$C2_Invasiveness, model_weights)
ref_lat <- weighted.mean(df_agg$wcg_lat_scaled, model_weights)
ref_lng <- weighted.mean(df_agg$wcg_lng_scaled, model_weights)
ref_form <- levels(df_agg$propagule_type)[1]
ref_year <- levels(df_agg$year_label)[1]

my_theme <- theme_classic() +
  theme(
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = NA),
    panel.grid.major = element_line(colour = "white", linetype = "solid"),
    text = element_text(family = "serif"),
    plot.title = element_text(family = "serif", colour = "black", size = 14),
    plot.subtitle = element_text(family = "serif", colour = "black", size = 14),
    axis.title = element_text(family = "serif", colour = "black", size = 14),
    axis.text = element_text(family = "serif", colour = "black", size = 14),
    legend.position = c(1.0, 0.01),
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.title = element_text(family = "serif", size = 14),
    legend.text = element_text(family = "serif", size = 14)
  )
my_main_theme <- my_theme + theme(legend.position = "none")
# A. Naturalization
dat_m1 <- ggemmeans(
  best_model,
  terms = "C1_NatStatus [-0.6666667, 0.3333333]",
  condition = list(C2_Invasiveness = 0,
                   propagule_type = ref_form,
                   wcg_lat_scaled = ref_lat,
                   wcg_lng_scaled = ref_lng,
                   year_label = ref_year)
)
dat_m1$x <- factor(dat_m1$x,
                   labels = c("Non-naturalized", "Naturalized"))
p6 <- ggplot(dat_m1, aes(x = x, y = predicted, fill = x)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  scale_x_discrete() +
  scale_fill_aaas() +
  scale_y_continuous(labels = scales::label_number(scale = 100), limits = c(0, 1)) +
  labs(title = "Naturalization: ***", x = NULL, y = "Percent of taxa\noutside climatic ranges") +
  my_main_theme

# B. Invasiveness
dat_m2 <- ggemmeans(
  best_model,
  terms = "C2_Invasiveness [-0.5, 0.5]",
  condition = list(C1_NatStatus = 1/3,
                   propagule_type = ref_form,
                   wcg_lat_scaled = ref_lat,
                   wcg_lng_scaled = ref_lng,
                   year_label = ref_year)
)
p7 <- ggplot(dat_m2, aes(x = factor(x), y = predicted, fill = factor(x))) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  scale_x_discrete(labels = c("-0.5" = "Non-invasive", "0.5" = "Invasive")) +
  scale_fill_aaas() +
  scale_y_continuous(labels = scales::label_number(scale = 100), limits = c(0, 1)) +
  labs(title = "Invasiveness: ***", x = NULL, y = "Percent of taxa\noutside climatic ranges") +
  my_main_theme

# C. Sale form
dat_m3 <- ggemmeans(
  best_model,
  terms = "propagule_type",
  condition = list(C1_NatStatus = ref_nat,
                   C2_Invasiveness = ref_invasive,
                   wcg_lat_scaled = ref_lat,
                   wcg_lng_scaled = ref_lng,
                   year_label = ref_year)
)
dat_m3$x <- factor(dat_m3$x, levels = c("seed", "adult"))

p8 <- ggplot(dat_m3, aes(x = x, y = predicted, fill = x)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  scale_x_discrete(labels = c("seed" = "Seeds", "adult" = "Live plants")) +
  scale_fill_aaas() +
  scale_y_continuous(labels = scales::label_number(scale = 100), limits = c(0, 1)) +
  labs(title = "Sale form: ns", x = NULL, y = "Percent of taxa\noutside climatic ranges") +
  my_main_theme

# D. Year
dat_m4 <- ggemmeans(
  best_model,
  terms = "year_label",
  condition = list(C1_NatStatus = ref_nat,
                   C2_Invasiveness = ref_invasive,
                   propagule_type = ref_form,
                   wcg_lat_scaled = ref_lat,
                   wcg_lng_scaled = ref_lng)
)
p9 <- ggplot(dat_m4, aes(x = x, y = predicted,  fill = x)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  scale_fill_aaas() +
  scale_y_continuous(labels = scales::label_number(scale = 100), limits = c(0, 1)) +
  labs(title = "Year: ns", x = NULL, y = "Percent of taxa\noutside climatic ranges") +
  my_main_theme

aaas_blue <- pal_aaas("default")(10)[1]
# E. Latitude
dat_m5 <- ggemmeans(
  best_model,
  terms = "wcg_lat_scaled",
  condition = list(C1_NatStatus = ref_nat,
                   C2_Invasiveness = ref_invasive,
                   propagule_type = ref_form,
                   wcg_lng_scaled = ref_lng,
                   year_label = ref_year)
)
p10 <- plot(dat_m5) +
  geom_line(color = aaas_blue, linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              fill = aaas_blue, alpha = 0.15) +
  labs(title = "Latitude: ***",
       x = "Standardized latitude",
       y = "Percent of taxa\noutside climatic ranges") +
  scale_y_continuous(labels = scales::label_number(scale = 100), limits = c(0, 1)) +
  my_main_theme

# F. Longitude
dat_m6 <- ggemmeans(
  best_model,
  terms = "wcg_lng_scaled",
  condition = list(C1_NatStatus = ref_nat,
                   C2_Invasiveness = ref_invasive,
                   propagule_type = ref_form,
                   wcg_lat_scaled = ref_lat,
                   year_label = ref_year)
)
p11 <- plot(dat_m6) +
  geom_line(color = aaas_blue, linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              fill = aaas_blue, alpha = 0.15) +
  labs(title = "Longitude: ***",
       x = "Standardized longitude",
       y = "Percent of taxa\noutside climatic ranges") +
  scale_y_continuous(labels = scales::label_number(scale = 100), limits = c(0, 1)) +
  my_main_theme

plot_mian <- wrap_plots(p6, p7, p8, p9, p10, p11, nrow = 3) +
  plot_annotation(tag_levels = 'a')

ggsave("results/Figure_2.png", plot_mian, device = ragg::agg_png,
       width = 8, height = 10, units = "in", dpi = 600, bg = "white")
ggsave("results/Figure_2.pdf", plot_mian, device = grDevices::cairo_pdf,
       width = 8, height = 10, units = "in", family = "serif", bg = "white")
