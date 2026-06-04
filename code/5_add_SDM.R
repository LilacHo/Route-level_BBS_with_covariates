# 5_add_SDM.R
# Add rcp45 and rcp85 columns to output/CH_no_habitat_grassland_routes.csv
# by extracting SDM classified-change raster values at each route location,
# per species.
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

# Read the route-level CH_no_habitat output --------------------------------
ch_routes <- read.csv(here::here("output", "CH_no_habitat_grassland_routes.csv"))

cat("Input rows:", nrow(ch_routes), "\n")
cat("Species:", length(unique(ch_routes$species_code)), "\n")

# Initialize new columns
ch_routes$rcp45 <- NA
ch_routes$rcp85 <- NA

# Loop through each species and extract raster values ----------------------
species_codes <- unique(ch_routes$species_code)

for (abbr in species_codes) {
  cat("Processing:", abbr, "\n")

  # Subset routes for this species
  idx <- which(ch_routes$species_code == abbr)
  sp_routes <- ch_routes[idx, ]

  # Build raster file paths
  rcp45_file <- here::here("data", "rcp45_grasslands", abbr,
                           paste0("grasslands_", abbr, "_breeding_2025_45_ENSEMBLE_classifiedchange.tif"))
  rcp85_file <- here::here("data", "rcp85_grasslands", abbr,
                           paste0("grasslands_", abbr, "_breeding_2025_85_ENSEMBLE_classifiedchange.tif"))

  # Check that raster files exist

  if (!file.exists(rcp45_file)) {
    message("  WARNING: rcp45 raster not found for ", abbr, " — skipping")
    next
  }
  if (!file.exists(rcp85_file)) {
    message("  WARNING: rcp85 raster not found for ", abbr, " — skipping")
    next
  }

  # Load rasters
  rcp45_rast <- rast(rcp45_file)
  rcp85_rast <- rast(rcp85_file)

  # Convert route locations to SpatVector, project to raster CRS
  routes_sv <- vect(sp_routes, geom = c("longitude", "latitude"), crs = "EPSG:4326")
  routes_sv_proj <- project(routes_sv, crs(rcp45_rast))

  # Extract raster values at route points
  vals_45 <- extract(rcp45_rast, routes_sv_proj)
  vals_85 <- extract(rcp85_rast, routes_sv_proj)

  # Assign extracted values back (column 2 holds the raster values)
  ch_routes$rcp45[idx] <- vals_45[, 2]
  ch_routes$rcp85[idx] <- vals_85[, 2]
}

# Summary ------------------------------------------------------------------
cat("\n=== Summary ===\n")
cat("Total rows:", nrow(ch_routes), "\n")
cat("rcp45 non-NA:", sum(!is.na(ch_routes$rcp45)), "\n")
cat("rcp85 non-NA:", sum(!is.na(ch_routes$rcp85)), "\n\n")

print(ch_routes %>%
        group_by(species_code) %>%
        summarise(n_routes = n(),
                  rcp45_mean = mean(rcp45, na.rm = TRUE),
                  rcp85_mean = mean(rcp85, na.rm = TRUE),
                  .groups = "drop"))

# Overwrite the output CSV with new columns --------------------------------
write.csv(ch_routes,
          here::here("output", "CH_no_habitat_grassland_routes.csv"),
          row.names = FALSE)

cat("\nWrote: output/CH_no_habitat_grassland_routes.csv (with rcp45 and rcp85 columns)\n")
