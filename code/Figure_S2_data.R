# Figure S2: calculate richness and suitable areas from source rasters.
library(terra)
library(readr)

output_dir <- "results"
dir.create(output_dir, showWarnings = FALSE)
classes <- c("inv", "nat", "nonnat")
class_english <- c(inv = "Invasive", nat = "Naturalized non-invasive",
                   nonnat = "Non-naturalized")

# 1. Match taxa to rasters.
records <- as.data.frame(read_csv(
  "data/Step4_Refugia_Full_Filtered_NotVeg.csv",
  col_types = cols(.default = col_character()), show_col_types = FALSE),
  stringsAsFactors = FALSE)
records$species_dir <- trimws(records$species_dir)
records$invasion_label <- trimws(records$invasion_label)
species <- unique(records[c("species_dir", "invasion_label")])
species <- species[order(species$species_dir), , drop = FALSE]
rownames(species) <- NULL
stopifnot(!anyDuplicated(species$species_dir), nrow(species) == 748L,
          setequal(unique(species$invasion_label), classes))
record_counts <- table(records$species_dir)
species$csv_record_count <- as.integer(record_counts[species$species_dir])

expected_counts <- unique(records[c("species_dir", "presence_cells", "absence_cells")])
stopifnot(!anyDuplicated(expected_counts$species_dir))
expected_counts$presence_cells <- as.numeric(expected_counts$presence_cells)
expected_counts$absence_cells <- as.numeric(expected_counts$absence_cells)
count_index <- match(species$species_dir, expected_counts$species_dir)

files <- list.files("data/SDM_combined", pattern = "_combined_binary[.]tif$",
                    recursive = TRUE, full.names = TRUE)
files <- files[basename(dirname(files)) == "EnsemblePlot"]
identifiers <- sub("_combined_binary[.]tif$", "", basename(files))
stopifnot(!anyDuplicated(identifiers))
species$source_path <- files[match(species$species_dir, identifiers)]
stopifnot(!anyNA(species$source_path))
species$file_status <- rep("matched", nrow(species))
extra <- !identifiers %in% species$species_dir

# 2. Cell areas and richness.
template <- rast(species$source_path[1])
stopifnot(nlyr(template) == 1L, is.lonlat(template))
domain <- !is.na(values(template, mat = FALSE))
cell_area <- cellSize(template, mask = FALSE, unit = "km")
areas <- values(cell_area, mat = FALSE)
domain_area <- sum(areas[domain])
richness <- matrix(0L, nrow = ncell(template), ncol = length(classes),
                   dimnames = list(NULL, classes))
species$presence_cells <- species$absence_cells <- species$nodata_cells <- 0L
species$suitable_area_km2 <- species$valid_area_km2 <- NA_real_
species$mask_difference_cells <- 0L
species$source_size_bytes <- file.info(species$source_path)$size
species$source_md5 <- unname(tools::md5sum(species$source_path))
species$selected_copy <- NA_character_
species$copy_md5_match <- NA
species$csv_presence_cells <- expected_counts$presence_cells[count_index]
species$csv_absence_cells <- expected_counts$absence_cells[count_index]

for (i in seq_len(nrow(species))) {
  raster <- rast(species$source_path[i])
  stopifnot(nlyr(raster) == 1L, compareGeom(template, raster))
  raster_values <- values(raster, mat = FALSE)
  valid <- !is.na(raster_values)
  stopifnot(all(raster_values[valid] %in% c(0, 1)))
  present <- valid & raster_values == 1
  class_index <- match(species$invasion_label[i], classes)
  richness[, class_index] <- richness[, class_index] + as.integer(present)
  species$presence_cells[i] <- sum(present)
  species$absence_cells[i] <- sum(valid & raster_values == 0)
  species$nodata_cells[i] <- sum(!valid)
  species$suitable_area_km2[i] <- sum(areas[present])
  species$valid_area_km2[i] <- sum(areas[valid])
  species$mask_difference_cells[i] <- sum(valid != domain)
}
species$csv_counts_match <- species$presence_cells == species$csv_presence_cells &
  species$absence_cells == species$csv_absence_cells
stopifnot(all(species$mask_difference_cells == 0L))
richness[!domain, ] <- NA_integer_
total_richness <- rowSums(richness)
class_presence <- richness > 0
overlap_code <- as.integer(class_presence[, 1]) +
  2L * as.integer(class_presence[, 2]) + 4L * as.integer(class_presence[, 3])

# 3. Summarize areas and group overlap.
# SD measures variation among taxa.
summary_rows <- vector("list", length(classes))
for (j in seq_along(classes)) {
  label <- classes[j]
  present <- !is.na(overlap_code) & class_presence[, j]
  exclusive <- !is.na(overlap_code) & overlap_code == 2^(j - 1L)
  species_areas <- species$suitable_area_km2[species$invasion_label == label]
  summary_rows[[j]] <- data.frame(
    invasion_label = label,
    class_en = unname(class_english[label]), species_count = length(species_areas),
    occupied_cells = sum(present), union_area_km2 = sum(areas[present]),
    coverage_pct_of_valid_grid = 100 * sum(areas[present]) / domain_area,
    exclusive_area_km2 = sum(areas[exclusive]),
    max_species_per_grid = max(richness[, j], na.rm = TRUE),
    area_weighted_mean_species_per_grid = sum(richness[domain, j] * areas[domain]) / domain_area,
    summed_species_area_km2 = sum(species_areas),
    mean_species_area_km2 = mean(species_areas),
    median_species_area_km2 = stats::median(species_areas),
    sd_species_area_km2 = stats::sd(species_areas),
    min_species_area_km2 = min(species_areas),
    q25_species_area_km2 = unname(stats::quantile(species_areas, 0.25)),
    q75_species_area_km2 = unname(stats::quantile(species_areas, 0.75)),
    max_species_area_km2 = max(species_areas), zero_area_species = sum(species_areas == 0),
    stringsAsFactors = FALSE)
}
summary <- do.call(rbind, summary_rows)
rownames(summary) <- NULL
mean_columns <- c("invasion_label", "class_en", "species_count",
                  "mean_species_area_km2", "median_species_area_km2", "sd_species_area_km2",
                  "min_species_area_km2", "q25_species_area_km2", "q75_species_area_km2",
                  "max_species_area_km2", "zero_area_species")
mean_area <- summary[mean_columns]
mean_area$mean_pct_of_valid_grid <- 100 * mean_area$mean_species_area_km2 / domain_area

overlap <- data.frame(
  bitmask = 0:7, inv_present = as.integer(bitwAnd(0:7, 1L) != 0L),
  nat_present = as.integer(bitwAnd(0:7, 2L) != 0L),
  nonnat_present = as.integer(bitwAnd(0:7, 4L) != 0L),
  combination = c("none", "inv_only", "nat_only", "inv_nat", "nonnat_only",
                   "inv_nonnat", "nat_nonnat", "inv_nat_nonnat"),
  cells = 0L, area_km2 = 0)
for (i in seq_len(nrow(overlap))) {
  selected <- !is.na(overlap_code) & overlap_code == overlap$bitmask[i]
  overlap$cells[i] <- sum(selected)
  overlap$area_km2[i] <- sum(areas[selected])
}
overlap$pct_of_valid_grid <- 100 * overlap$area_km2 / domain_area

# 4. Area totals and grid data.
stopifnot(sum(overlap$cells) == sum(domain))
stopifnot(abs(sum(overlap$area_km2) - domain_area) < 1e-6)
for (j in seq_along(classes)) {
  selected <- overlap[[paste0(classes[j], "_present")]] == 1L
  stopifnot(abs(sum(overlap$area_km2[selected]) - summary$union_area_km2[j]) < 1e-6)
  stopifnot(abs(sum(richness[domain, j] * areas[domain]) -
                summary$summed_species_area_km2[j]) < 1e-4)
}
stopifnot(all(abs(summary$mean_species_area_km2 * summary$species_count -
                   summary$summed_species_area_km2) < 1e-4))

grid <- as.data.frame(terra::xyFromCell(template, which(domain)))
grid$grid_id <- which(domain)
grid$cell_area_km2 <- areas[domain]
for (label in classes) grid[[paste0(label, "_species_count")]] <- richness[domain, label]
grid$all_species_count <- total_richness[domain]
grid$overlap_bitmask <- overlap_code[domain]
overall <- data.frame(
  csv_rows = nrow(records), unique_species = nrow(species),
  selected_combined_files = nrow(species), available_combined_files = length(files),
  extra_files_excluded = sum(extra), excluded_species = paste(sort(identifiers[extra]), collapse = "; "),
  valid_grid_cells = sum(domain), valid_grid_area_km2 = domain_area,
  all_union_cells = sum(total_richness > 0, na.rm = TRUE),
  all_union_area_km2 = sum(areas[!is.na(total_richness) & total_richness > 0]),
  all_max_species_per_grid = max(total_richness, na.rm = TRUE),
  csv_pixel_count_mismatches = sum(!species$csv_counts_match),
  differing_masks = sum(species$mask_difference_cells > 0),
  all_copies_md5_verified = NA,
  nrow = terra::nrow(template), ncol = terra::ncol(template),
  x_resolution_degrees = terra::res(template)[1], y_resolution_degrees = terra::res(template)[2],
  crs = terra::crs(template, proj = TRUE), csv_md5 = unname(tools::md5sum("data/Step4_Refugia_Full_Filtered_NotVeg.csv")),
  r_version = R.version.string, terra_version = as.character(utils::packageVersion("terra")),
  readr_version = as.character(utils::packageVersion("readr")), stringsAsFactors = FALSE)

# 5. Save rasters and tables.
raster_output <- file.path(output_dir, "rasters")
dir.create(raster_output, showWarnings = FALSE)
stopifnot(nrow(species) < 65535L)
for (label in classes) {
  output_name <- paste0(label, "_species_richness")
  output_raster <- setValues(rast(template), richness[, label])
  names(output_raster) <- output_name
  writeRaster(output_raster, file.path(raster_output, paste0(output_name, ".tif")),
               overwrite = TRUE,
               wopt = list(datatype = "INT2U", NAflag = 65535, gdal = "COMPRESS=DEFLATE"))
}
results <- list(grid = grid, summary = summary, mean_area = mean_area, overall = overall,
                extent = as.vector(ext(template)), resolution = res(template), overlap = overlap)
write_excel_csv(species, "results/species_areas_and_checks.csv")
write_excel_csv(mean_area, "results/class_mean_species_area.csv")
write_excel_csv(grid, "results/grid_species_counts.csv")
write_excel_csv(overall, "results/overall_summary_and_qa.csv")
saveRDS(results, "data/plot_data.rds")
