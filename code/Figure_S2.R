# Figure S2: taxon richness and mean suitable area.
library(ggplot2)
library(patchwork)
library(sf)

dir.create("results", showWarnings = FALSE)
results <- readRDS("data/plot_data.rds")
grid_data <- results$grid
area_summary <- results$summary
font_scale <- 16 / 10

# 1. Map boundaries.
china_raw <- st_read("data/China_map/bou1_4p.shp", options = "ENCODING=GBK", quiet = TRUE)
south <- st_read("data/China_map/south_sea.shp", quiet = TRUE)
st_crs(china_raw) <- 4326
china_raw <- st_transform(china_raw, 4326)
south <- st_transform(south, 4326)
dash <- south[south$GBCODE == 61010 & south$BOU1_4M_ != 691, ]
stopifnot(nrow(dash) == 9L)
islets <- south[south$GBCODE != 61010, ]
land_geometry <- st_union(st_make_valid(st_set_crs(st_geometry(china_raw), NA)))
frame <- st_as_sfc(st_bbox(c(xmin = 70, ymin = 1, xmax = 145, ymax = 57)))
outside_geometry <- st_difference(frame, land_geometry)
land <- st_sf(geometry = st_set_crs(land_geometry, 4326))
outside <- st_sf(geometry = st_set_crs(outside_geometry, 4326))

land_path <- as.data.frame(st_coordinates(land))
land_path$group <- paste(1L, land_path$L2, land_path$L3, sep = "_")
land_path$subgroup <- paste(1L, land_path$L1, land_path$L2, land_path$L3, sep = "_")
outside_path <- as.data.frame(st_coordinates(outside))
outside_path$group <- paste(1L, outside_path$L2, sep = "_")
outside_path$subgroup <- paste(1L, outside_path$L1, outside_path$L2, sep = "_")
islet_path <- as.data.frame(st_coordinates(islets))
islet_path$group <- paste(1L, islet_path$L1, sep = "_")
islet_path$subgroup <- islet_path$group
dash_path <- as.data.frame(st_coordinates(dash))
dash_path$group <- paste(1L, dash_path$L1, sep = "_")
dash_path$subgroup <- dash_path$group

# 2. Colours, text and raster grid.
classes <- c("nonnat", "nat", "inv")
labels <- c(inv = "Invasive", nat = "Naturalized non-invasive", nonnat = "Non-naturalized")
bar_colors <- c(inv = "#547EA0", nat = "#547EA0", nonnat = "#547EA0")
continuous <- c("#E0F1DA", "#B7DFAB", "#78C679", "#31A354", "#08783E", "#004D2A")
zero_grey <- "#C9C9C9"
nodata_grey <- "#E9E9E9"
publication_theme <- theme_classic(base_size = 16, base_family = "serif") +
  theme(axis.title = element_text(size = 10 * font_scale),
    axis.text = element_text(size = 9 * font_scale, colour = "black"),
    axis.line = element_line(linewidth = .3),
    axis.ticks = element_line(linewidth = .3),
    axis.ticks.length = grid::unit(1.2, "mm"),
    plot.title = element_text(size = 10.5 * font_scale, face = "plain", hjust = .5),
    plot.subtitle = element_text(size = 9.5 * font_scale, hjust = .5),
    plot.tag = element_text(size = 16, face = "plain"),
    plot.tag.position = "topleft",
    legend.title = element_text(size = 10 * font_scale),
    legend.text = element_text(size = 9 * font_scale),
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(2, 2, 2, 2, "mm"))
number_of_rows <- results$overall$nrow
number_of_columns <- results$overall$ncol
full_grid <- expand.grid(
  x = results$extent[1] + (seq_len(number_of_columns) - .5) * results$resolution[1],
  y = results$extent[4] - (seq_len(number_of_rows) - .5) * results$resolution[2])
common_max <- max(grid_data[, paste0(classes, "_species_count")])


# A. Non-naturalized.
raster_nonnat <- full_grid
raster_nonnat$richness <- NA_real_
raster_nonnat$richness[grid_data$grid_id] <- grid_data$nonnat_species_count
raster_nonnat$is_zero <- as.numeric(!is.na(raster_nonnat$richness) & raster_nonnat$richness == 0)
raster_nonnat$positive_richness <- raster_nonnat$richness
raster_nonnat$positive_richness[which(raster_nonnat$richness == 0)] <- NA_real_

map_nonnat <- ggplot() +
  geom_polygon(data = land_path, aes(X, Y, group = group, subgroup = subgroup),
    rule = "evenodd", fill = nodata_grey, colour = NA) +
  geom_raster(data = raster_nonnat, aes(x, y, alpha = is_zero), fill = zero_grey,
    interpolate = FALSE, show.legend = FALSE) +
  scale_alpha_identity(guide = "none") +
  geom_raster(data = raster_nonnat, aes(x, y, fill = positive_richness),
    interpolate = FALSE, show.legend = TRUE) +
  scale_fill_gradientn(colours = continuous, limits = c(1, common_max),
    breaks = c(1, 50, 100, 150, 200, 250, common_max),
    name = "Predicted taxon richness", na.value = "transparent",
    guide = guide_colourbar(title.position = "top",
      title.hjust = .5, barwidth = grid::unit(125, "mm"),
      barheight = grid::unit(3.5, "mm"), ticks.colour = "white")) +
  geom_polygon(data = outside_path, aes(X, Y, group = group, subgroup = subgroup),
    rule = "evenodd", fill = "white", colour = NA) +
  geom_polygon(data = land_path, aes(X, Y, group = group, subgroup = subgroup),
    rule = "evenodd", fill = NA, colour = "#4D5658", linewidth = .28) +
  geom_path(data = islet_path, aes(X, Y, group = group),
    colour = "#586260", linewidth = .15) +
  geom_path(data = dash_path, aes(X, Y, group = group),
    colour = "#313A3B", linewidth = .35) +
  labs(x = NULL, y = NULL) +
  coord_fixed(ratio = 1 / cos(35.5 * pi / 180),
    xlim = c(72, 143), ylim = c(16, 55), expand = FALSE) +
  scale_x_continuous(breaks = c(80, 100, 120), labels = paste0(c(80, 100, 120), "\u00B0E")) +
  scale_y_continuous(breaks = c(20, 35, 50), labels = paste0(c(20, 35, 50), "\u00B0N")) +
  publication_theme +
  theme(legend.position = "bottom", axis.line = element_blank(),
    panel.border = element_rect(fill = NA, colour = "#AEB5B5", linewidth = .3))
map_nonnat <- map_nonnat + labs(tag = "a", title = labels["nonnat"],
  subtitle = paste0("n = ", area_summary$species_count[area_summary$invasion_label == "nonnat"]))

inset_nonnat <- ggplot() +
  geom_polygon(data = land_path, aes(X, Y, group = group, subgroup = subgroup),
    rule = "evenodd", fill = nodata_grey, colour = NA) +
  geom_raster(data = raster_nonnat, aes(x, y, alpha = is_zero), fill = zero_grey,
    interpolate = FALSE, show.legend = FALSE) +
  scale_alpha_identity(guide = "none") +
  geom_raster(data = raster_nonnat, aes(x, y, fill = positive_richness),
    interpolate = FALSE, show.legend = FALSE) +
  scale_fill_gradientn(colours = continuous, limits = c(1, common_max),
    breaks = c(1, 50, 100, 150, 200, 250, common_max),
    name = "Predicted taxon richness", na.value = "transparent",
    guide = "none") +
  geom_polygon(data = outside_path, aes(X, Y, group = group, subgroup = subgroup),
    rule = "evenodd", fill = "white", colour = NA) +
  geom_polygon(data = land_path, aes(X, Y, group = group, subgroup = subgroup),
    rule = "evenodd", fill = NA, colour = "#4D5658", linewidth = .23) +
  geom_path(data = islet_path, aes(X, Y, group = group),
    colour = "#586260", linewidth = .23) +
  geom_path(data = dash_path, aes(X, Y, group = group),
    colour = "#313A3B", linewidth = .55) +
  labs(x = NULL, y = NULL) +
  coord_fixed(ratio = 1 / cos(13.5 * pi / 180),
    xlim = c(105, 125.5), ylim = c(3, 24), expand = FALSE) +
  theme_void(base_family = "serif") +
  theme(panel.border = element_rect(fill = NA, colour = "#7D8484", linewidth = .3),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(0, 0, 0, 0), legend.position = "none")
panel_nonnat <- map_nonnat + inset_element(inset_nonnat,
  left = .79, bottom = .025, right = .995, top = .40,
  align_to = "panel", on_top = TRUE)

# B. Naturalized non-invasive.
raster_nat <- full_grid
raster_nat$richness <- NA_real_
raster_nat$richness[grid_data$grid_id] <- grid_data$nat_species_count
raster_nat$is_zero <- as.numeric(!is.na(raster_nat$richness) & raster_nat$richness == 0)
raster_nat$positive_richness <- raster_nat$richness
raster_nat$positive_richness[which(raster_nat$richness == 0)] <- NA_real_

map_nat <- ggplot() +
  geom_polygon(data = land_path, aes(X, Y, group = group, subgroup = subgroup),
    rule = "evenodd", fill = nodata_grey, colour = NA) +
  geom_raster(data = raster_nat, aes(x, y, alpha = is_zero), fill = zero_grey,
    interpolate = FALSE, show.legend = FALSE) +
  scale_alpha_identity(guide = "none") +
  geom_raster(data = raster_nat, aes(x, y, fill = positive_richness),
    interpolate = FALSE, show.legend = TRUE) +
  scale_fill_gradientn(colours = continuous, limits = c(1, common_max),
    breaks = c(1, 50, 100, 150, 200, 250, common_max),
    name = "Predicted taxon richness", na.value = "transparent",
    guide = guide_colourbar(title.position = "top",
      title.hjust = .5, barwidth = grid::unit(125, "mm"),
      barheight = grid::unit(3.5, "mm"), ticks.colour = "white")) +
  geom_polygon(data = outside_path, aes(X, Y, group = group, subgroup = subgroup),
    rule = "evenodd", fill = "white", colour = NA) +
  geom_polygon(data = land_path, aes(X, Y, group = group, subgroup = subgroup),
    rule = "evenodd", fill = NA, colour = "#4D5658", linewidth = .28) +
  geom_path(data = islet_path, aes(X, Y, group = group),
    colour = "#586260", linewidth = .15) +
  geom_path(data = dash_path, aes(X, Y, group = group),
    colour = "#313A3B", linewidth = .35) +
  labs(x = NULL, y = NULL) +
  coord_fixed(ratio = 1 / cos(35.5 * pi / 180),
    xlim = c(72, 143), ylim = c(16, 55), expand = FALSE) +
  scale_x_continuous(breaks = c(80, 100, 120), labels = paste0(c(80, 100, 120), "\u00B0E")) +
  scale_y_continuous(breaks = c(20, 35, 50), labels = paste0(c(20, 35, 50), "\u00B0N")) +
  publication_theme +
  theme(legend.position = "bottom", axis.line = element_blank(),
    panel.border = element_rect(fill = NA, colour = "#AEB5B5", linewidth = .3))
map_nat <- map_nat + labs(tag = "b", title = labels["nat"],
  subtitle = paste0("n = ", area_summary$species_count[area_summary$invasion_label == "nat"]))

inset_nat <- ggplot() +
  geom_polygon(data = land_path, aes(X, Y, group = group, subgroup = subgroup),
    rule = "evenodd", fill = nodata_grey, colour = NA) +
  geom_raster(data = raster_nat, aes(x, y, alpha = is_zero), fill = zero_grey,
    interpolate = FALSE, show.legend = FALSE) +
  scale_alpha_identity(guide = "none") +
  geom_raster(data = raster_nat, aes(x, y, fill = positive_richness),
    interpolate = FALSE, show.legend = FALSE) +
  scale_fill_gradientn(colours = continuous, limits = c(1, common_max),
    breaks = c(1, 50, 100, 150, 200, 250, common_max),
    name = "Predicted taxon richness", na.value = "transparent",
    guide = "none") +
  geom_polygon(data = outside_path, aes(X, Y, group = group, subgroup = subgroup),
    rule = "evenodd", fill = "white", colour = NA) +
  geom_polygon(data = land_path, aes(X, Y, group = group, subgroup = subgroup),
    rule = "evenodd", fill = NA, colour = "#4D5658", linewidth = .23) +
  geom_path(data = islet_path, aes(X, Y, group = group),
    colour = "#586260", linewidth = .23) +
  geom_path(data = dash_path, aes(X, Y, group = group),
    colour = "#313A3B", linewidth = .55) +
  labs(x = NULL, y = NULL) +
  coord_fixed(ratio = 1 / cos(13.5 * pi / 180),
    xlim = c(105, 125.5), ylim = c(3, 24), expand = FALSE) +
  theme_void(base_family = "serif") +
  theme(panel.border = element_rect(fill = NA, colour = "#7D8484", linewidth = .3),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(0, 0, 0, 0), legend.position = "none")
panel_nat <- map_nat + inset_element(inset_nat,
  left = .79, bottom = .025, right = .995, top = .40,
  align_to = "panel", on_top = TRUE)

# C. Invasive.
raster_inv <- full_grid
raster_inv$richness <- NA_real_
raster_inv$richness[grid_data$grid_id] <- grid_data$inv_species_count
raster_inv$is_zero <- as.numeric(!is.na(raster_inv$richness) & raster_inv$richness == 0)
raster_inv$positive_richness <- raster_inv$richness
raster_inv$positive_richness[which(raster_inv$richness == 0)] <- NA_real_

map_inv <- ggplot() +
  geom_polygon(data = land_path, aes(X, Y, group = group, subgroup = subgroup),
    rule = "evenodd", fill = nodata_grey, colour = NA) +
  geom_raster(data = raster_inv, aes(x, y, alpha = is_zero), fill = zero_grey,
    interpolate = FALSE, show.legend = FALSE) +
  scale_alpha_identity(guide = "none") +
  geom_raster(data = raster_inv, aes(x, y, fill = positive_richness),
    interpolate = FALSE, show.legend = TRUE) +
  scale_fill_gradientn(colours = continuous, limits = c(1, common_max),
    breaks = c(1, 50, 100, 150, 200, 250, common_max),
    name = "Predicted taxon richness", na.value = "transparent",
    guide = guide_colourbar(title.position = "top",
      title.hjust = .5, barwidth = grid::unit(125, "mm"),
      barheight = grid::unit(3.5, "mm"), ticks.colour = "white")) +
  geom_polygon(data = outside_path, aes(X, Y, group = group, subgroup = subgroup),
    rule = "evenodd", fill = "white", colour = NA) +
  geom_polygon(data = land_path, aes(X, Y, group = group, subgroup = subgroup),
    rule = "evenodd", fill = NA, colour = "#4D5658", linewidth = .28) +
  geom_path(data = islet_path, aes(X, Y, group = group),
    colour = "#586260", linewidth = .15) +
  geom_path(data = dash_path, aes(X, Y, group = group),
    colour = "#313A3B", linewidth = .35) +
  labs(x = NULL, y = NULL) +
  coord_fixed(ratio = 1 / cos(35.5 * pi / 180),
    xlim = c(72, 143), ylim = c(16, 55), expand = FALSE) +
  scale_x_continuous(breaks = c(80, 100, 120), labels = paste0(c(80, 100, 120), "\u00B0E")) +
  scale_y_continuous(breaks = c(20, 35, 50), labels = paste0(c(20, 35, 50), "\u00B0N")) +
  publication_theme +
  theme(legend.position = "bottom", axis.line = element_blank(),
    panel.border = element_rect(fill = NA, colour = "#AEB5B5", linewidth = .3))
map_inv <- map_inv + labs(tag = "c", title = labels["inv"],
  subtitle = paste0("n = ", area_summary$species_count[area_summary$invasion_label == "inv"]))

inset_inv <- ggplot() +
  geom_polygon(data = land_path, aes(X, Y, group = group, subgroup = subgroup),
    rule = "evenodd", fill = nodata_grey, colour = NA) +
  geom_raster(data = raster_inv, aes(x, y, alpha = is_zero), fill = zero_grey,
    interpolate = FALSE, show.legend = FALSE) +
  scale_alpha_identity(guide = "none") +
  geom_raster(data = raster_inv, aes(x, y, fill = positive_richness),
    interpolate = FALSE, show.legend = FALSE) +
  scale_fill_gradientn(colours = continuous, limits = c(1, common_max),
    breaks = c(1, 50, 100, 150, 200, 250, common_max),
    name = "Predicted taxon richness", na.value = "transparent",
    guide = "none") +
  geom_polygon(data = outside_path, aes(X, Y, group = group, subgroup = subgroup),
    rule = "evenodd", fill = "white", colour = NA) +
  geom_polygon(data = land_path, aes(X, Y, group = group, subgroup = subgroup),
    rule = "evenodd", fill = NA, colour = "#4D5658", linewidth = .23) +
  geom_path(data = islet_path, aes(X, Y, group = group),
    colour = "#586260", linewidth = .23) +
  geom_path(data = dash_path, aes(X, Y, group = group),
    colour = "#313A3B", linewidth = .55) +
  labs(x = NULL, y = NULL) +
  coord_fixed(ratio = 1 / cos(13.5 * pi / 180),
    xlim = c(105, 125.5), ylim = c(3, 24), expand = FALSE) +
  theme_void(base_family = "serif") +
  theme(panel.border = element_rect(fill = NA, colour = "#7D8484", linewidth = .3),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(0, 0, 0, 0), legend.position = "none")
panel_inv <- map_inv + inset_element(inset_inv,
  left = .79, bottom = .025, right = .995, top = .40,
  align_to = "panel", on_top = TRUE)

# D. Mean suitable area.
area_summary$invasion_label <- factor(area_summary$invasion_label, levels = classes)
area_labels <- c(inv = "Invasive", nat = "Naturalized\nnon-invasive",
                 nonnat = "Non-\nnaturalized")
mean_area <- ggplot(area_summary, aes(x = invasion_label,
                           y = mean_species_area_km2 / 1e6, fill = invasion_label)) +
  geom_col(width = .58, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.2f", mean_species_area_km2 / 1e6)),
            vjust = -.6, size = 3.4 * font_scale, family = "serif") +
  scale_fill_manual(values = bar_colors, guide = "none") +
  scale_x_discrete(labels = area_labels) +
  scale_y_continuous(limits = c(0, 2.06), breaks = seq(0, 2, .5), expand = c(0, 0)) +
  labs(tag = "d", x = NULL,
       y = "Mean suitable area per taxon\n(10\u2076 km\u00B2)") +
  publication_theme +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 9 * font_scale, lineheight = 1.05),
        plot.margin = margin(2, 4, 2, 4, "mm"))
# 3. Arrange panels and export.
figure <- wrap_plots(panel_nonnat, panel_nat, panel_inv,
                    wrap_elements(full = mean_area), guide_area(),
                    design = "EE\nAB\nCD", heights = c(.34, 1, 1), guides = "collect") &
  theme(legend.position = "top")
# 8. Export PNG (600 dpi) and PDF.
dir.create(dirname("results/Figure_S2.png"), recursive = TRUE, showWarnings = FALSE)
ragg::agg_png("results/Figure_S2.png",
  width = 240, height = 205,
  units = "mm", res = 600, background = "white")
print(figure)
dev.off()
grDevices::cairo_pdf(file.path("results", "Figure_S2.pdf"),
  width = 240 / 25.4, height = 205 / 25.4,
  family = "serif", bg = "white")
print(figure)
invisible(dev.off())
