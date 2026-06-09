# 5_add_SDM.R
# Add rcp45 and rcp85 columns to each individual species CSV in
# output/species_routes/ by extracting SDM classified-change raster values
# at each route location. Writes revised CSVs to output/species_routes_sdm/.
#
# Species abbreviations are looked up from data/spp_names_codes_group_aou.csv
# using the species_code column already present in each route CSV.
#
# Raster value legend (from Smith et al. 2024):
#   0 = never suitable
#   1 = extirpation
#   2 = worsening (-50 to -100% change)
#   3 = slightly worsening (-25 to -50% change)
#   4 = neutral (-25 to +25% change)
#   5 = slightly improving (25 to 50% change)
#   6 = improving (50 to Inf % change)
#   7 = colonization

library(here)
library(tidyverse)
library(terra)
library(sf)

here::i_am("code/5_add_SDM.R")

# Read species name/code lookup table --------------------------------------
spp_df <- read.csv(here::here("data", "spp_names_codes_group_aou.csv"))

# List all species route CSVs ----------------------------------------------
species_routes_dir <- here::here("output", "species_routes")
species_files <- list.files(species_routes_dir, pattern = "_route_CH\\.csv$",
                            full.names = TRUE)

cat("Found", length(species_files), "species route files\n\n")

# Create output directory --------------------------------------------------
out_dir <- here::here("output", "species_routes_sdm")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Process each species file ------------------------------------------------
for (f in species_files) {

  # Skip if output already exists
  out_file <- file.path(out_dir, sub("_route_CH\\.csv$", "_route_CH_sdm.csv", basename(f)))
  if (file.exists(out_file)) {
    cat("Skipping (already exists):", basename(out_file), "\n")
    next
  }

  sp_routes <- read.csv(f)

  # Get the species abbreviation from the CSV (should be unique per file)
  abbr <- unique(sp_routes$species_code)
  if (length(abbr) != 1) {
    message("WARNING: Multiple species codes in ", basename(f), " — skipping")
    next
  }

  # Verify abbreviation exists in lookup table
  if (!abbr %in% spp_df$Code) {
    message("WARNING: Code '", abbr, "' not found in spp_names_codes_group_aou.csv — skipping")
    next
  }

  cat("Processing:", abbr, "(", basename(f), ")\n")

  # Initialize new columns

  sp_routes$rcp45 <- NA
  sp_routes$rcp85 <- NA

  # Build raster file paths
  rcp45_file <- here::here("data", "rcp45_grasslands", abbr,
                           paste0("grasslands_", abbr, "_breeding_2025_45_ENSEMBLE_classifiedchange.tif"))
  rcp85_file <- here::here("data", "rcp85_grasslands", abbr,
                           paste0("grasslands_", abbr, "_breeding_2025_85_ENSEMBLE_classifiedchange.tif"))

  # Check that raster files exist
  if (!file.exists(rcp45_file)) {
    message("  WARNING: rcp45 raster not found for ", abbr, " — skipping extraction")
  }
  if (!file.exists(rcp85_file)) {
    message("  WARNING: rcp85 raster not found for ", abbr, " — skipping extraction")
  }

  # Extract rcp45 values
  if (file.exists(rcp45_file)) {
    rcp45_rast <- rast(rcp45_file)
    routes_sv <- vect(sp_routes, geom = c("longitude", "latitude"), crs = "EPSG:4326")
    routes_sv_proj <- project(routes_sv, crs(rcp45_rast))
    vals_45 <- extract(rcp45_rast, routes_sv_proj)
    sp_routes$rcp45 <- vals_45[, 2]
  }

  # Extract rcp85 values
  if (file.exists(rcp85_file)) {
    rcp85_rast <- rast(rcp85_file)
    routes_sv <- vect(sp_routes, geom = c("longitude", "latitude"), crs = "EPSG:4326")
    routes_sv_proj <- project(routes_sv, crs(rcp85_rast))
    vals_85 <- extract(rcp85_rast, routes_sv_proj)
    sp_routes$rcp85 <- vals_85[, 2]
  }

  # Write revised CSV to output/species_routes_sdm/
  write.csv(sp_routes, out_file, row.names = FALSE)

  cat("  Wrote:", basename(out_file),
      "| rcp45 non-NA:", sum(!is.na(sp_routes$rcp45)),
      "| rcp85 non-NA:", sum(!is.na(sp_routes$rcp85)), "\n")
}

# Summary ------------------------------------------------------------------
cat("\n=== SDM extraction done ===\n")
cat("Output directory:", out_dir, "\n")
cat("Files written:", length(list.files(out_dir, pattern = "\\.csv$")), "\n")

# ==========================================================================
# Violin plots: mean CH_no_habitat by rcp45 category for each species
# ==========================================================================

# Raster category labels
rcp_labels <- c("0" = "Never suitable",
                "1" = "Extirpation",
                "2" = "Worsening",
                "3" = "Slightly worsening",
                "4" = "Neutral",
                "5" = "Slightly improving",
                "6" = "Improving",
                "7" = "Colonization")

# Create barplot output directory
plot_dir <- here::here("output", "species_routes_sdm_plot")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Read all SDM CSVs and generate bar plots
sdm_files <- list.files(out_dir, pattern = "_route_CH_sdm\\.csv$",
                        full.names = TRUE)

cat("\nGenerating bar plots for", length(sdm_files), "species...\n")

for (sf in sdm_files) {
  sp_data <- read.csv(sf)

  abbr <- unique(sp_data$species_code)
  sp_name <- unique(sp_data$species)

  # Skip if plot already exists
  plot_file <- file.path(plot_dir, paste0(abbr, "_CH_no_habitat_by_rcp.png"))
  if (file.exists(plot_file)) {
    cat("  Skipping plot (already exists):", basename(plot_file), "\n")
    next
  }

  # Skip if no rcp data at all
  if (all(is.na(sp_data$rcp45)) && all(is.na(sp_data$rcp85))) {
    cat("  Skipping plot for", abbr, "— no rcp values\n")
    next
  }

  # Reshape to long format for both scenarios
  long_45 <- sp_data %>%
    filter(!is.na(rcp45)) %>%
    mutate(category = factor(rcp45, levels = 0:7),
           scenario = "RCP 4.5") %>%
    select(category, scenario, CH_no_habitat)

  long_85 <- sp_data %>%
    filter(!is.na(rcp85)) %>%
    mutate(category = factor(rcp85, levels = 0:7),
           scenario = "RCP 8.5") %>%
    select(category, scenario, CH_no_habitat)

  plot_data <- bind_rows(long_45, long_85)

  # Create faceted violin plot (rcp45 and rcp85)
  p <- ggplot(plot_data, aes(x = category, y = CH_no_habitat, fill = scenario)) +
    geom_violin(trim = FALSE, scale = "width") +
    geom_boxplot(width = 0.1, outlier.size = 0.5, fill = "white", alpha = 0.6) +
    facet_wrap(~ scenario, ncol = 1) +
    scale_fill_manual(values = c("RCP 4.5" = "steelblue", "RCP 8.5" = "firebrick")) +
    labs(title = paste0(sp_name, " (", abbr, ")"),
         subtitle = paste0("n = ", nrow(sp_data), " routes"),
         x = "SDM classified change category",
         y = "Residual (%)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  # Save plot
  ggsave(plot_file, p, width = 8, height = 8, dpi = 150)
  cat("  Saved:", basename(plot_file), "\n")
}

cat("\n=== All done ===\n")
cat("Violin plots saved to:", plot_dir, "\n")
