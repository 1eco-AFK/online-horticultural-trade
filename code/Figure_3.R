# Figure 3: saved-model predictions.

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
# A. Latitude and naturalization
dat_c1_lat <- ggemmeans(
  best_model,
  terms = c("wcg_lat_scaled", "C1_NatStatus [-0.6666667, 0.3333333]"),
  condition = list(C2_Invasiveness = 0,
                   propagule_type = ref_form,
                   wcg_lng_scaled = ref_lng,
                   year_label = ref_year)
)
dat_c1_lat$group <- factor(dat_c1_lat$group,
                           labels = c("Non-naturalized", "Naturalized"))

p1 <- plot(dat_c1_lat,
           line_size = 1) +
  labs(title = "Latitude \u00D7 Naturalization: ***",
       y = "Percent of taxa\noutside climatic ranges", x = "Standardized latitude",
       color = "Naturalization status", fill = "Naturalization status") +
  scale_color_aaas() +
  scale_fill_aaas() +
  scale_y_continuous(
    labels = scales::label_number(scale = 100),
    limits = c(0, 1)
  ) +
  my_theme

# B. Longitude and naturalization
dat_c1_lng <- ggemmeans(
  best_model,
  terms = c("wcg_lng_scaled", "C1_NatStatus [-0.6666667, 0.3333333]"),
  condition = list(C2_Invasiveness = 0,
                   propagule_type = ref_form,
                   wcg_lat_scaled = ref_lat,
                   year_label = ref_year)
)
dat_c1_lng$group <- factor(dat_c1_lng$group,
                           labels = c("Non-naturalized", "Naturalized"))

p2 <- plot(dat_c1_lng,
           line_size = 1) +
  labs(title = "Longitude \u00D7 Naturalization: ***",
       y = "Percent of taxa\noutside climatic ranges", x = "Standardized longitude",
       color = "Naturalization status", fill = "Naturalization status") +
  scale_color_aaas() +
  scale_fill_aaas() +
  scale_y_continuous(
    labels = scales::label_number(scale = 100),
    limits = c(0, 1)
  ) +
  my_theme

# C. Latitude and invasiveness
dat_c2_lat <- ggemmeans(
  best_model,
  terms = c("wcg_lat_scaled", "C2_Invasiveness [-0.5, 0.5]"),
  condition = list(C1_NatStatus = 1/3,
                   propagule_type = ref_form,
                   wcg_lng_scaled = ref_lng,
                   year_label = ref_year)
)
dat_c2_lat$group <- as.factor(dat_c2_lat$group)

p3 <- plot(dat_c2_lat, line_size = 1) +
  labs(title = "Latitude \u00D7 Invasiveness: ***",
       y = "Percent of taxa\noutside climatic ranges", x = "Standardized latitude",
       color = "Invasiveness status", fill = "Invasiveness status") +
  scale_color_aaas(labels = c("-0.5" = "Non-invasive", "0.5" = "Invasive")) +
  scale_fill_aaas(labels = c("-0.5" = "Non-invasive", "0.5" = "Invasive")) +
  scale_y_continuous(labels = scales::label_number(scale = 100), limits = c(0, 1)) +
  my_theme

# D. Latitude and sale form
dat_prop_lat <- ggemmeans(
  best_model,
  terms = c("wcg_lat_scaled", "propagule_type"),
  condition = list(C1_NatStatus = ref_nat,
                   C2_Invasiveness = ref_invasive,
                   wcg_lng_scaled = ref_lng,
                   year_label = ref_year)
)
dat_prop_lat$group <- factor(dat_prop_lat$group,
                             levels = c("seed", "adult"),
                             labels = c("Seeds", "Live plants"))

p4 <- plot(dat_prop_lat, line_size = 1) +
  labs(title = "Latitude \u00D7 Sale form: ***",
       y = "Percent of taxa\noutside climatic ranges", x = "Standardized latitude",
       color = "Sale form", fill = "Sale form") +
  scale_color_aaas() +
  scale_fill_aaas() +
  scale_y_continuous(labels = scales::label_number(scale = 100), limits = c(0, 1)) +
  my_theme

plot_interaction <- wrap_plots(p1, p2, p3, p4, nrow = 2) +
  plot_annotation(tag_levels = 'A')

ggsave("results/Figure_3.png", plot_interaction, device = ragg::agg_png,
       width = 8, height = 8, units = "in", dpi = 600, bg = "white")
ggsave("results/Figure_3.pdf", plot_interaction, device = grDevices::cairo_pdf,
       width = 8, height = 8, units = "in", family = "serif", bg = "white")
