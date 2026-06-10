# 5_add_SDM.R
# Add rcp45 and rcp85 columns to each individual species CSV in
# output/species_routes/ by extracting SDM classified-change raster values
# at each route location. Writes revised CSVs to output/species_routes_sdm/.
#
# Species abbreviations are looked up from data/spp_names_codes_group_aou.csv
# using the species_code column already present in each route CSV.
#
# Raster value legend (from Bateman et al. 2020):
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

# ==========================================================================
# Combined: merge all species SDM files and create an overall violin plot
# ==========================================================================

cat("\nCombining all species SDM files...\n")

all_sdm <- sdm_files %>%
  map_dfr(read.csv)

cat("Total rows:", nrow(all_sdm), "\n")
cat("Species:", length(unique(all_sdm$species_code)), "\n")

# Overall violin plot — RCP 4.5 and RCP 8.5 faceted
long_all_45 <- all_sdm %>%
  filter(!is.na(rcp45)) %>%
  mutate(category = factor(rcp45, levels = 0:7),
         scenario = "RCP 4.5") %>%
  select(category, scenario, CH_no_habitat)

long_all_85 <- all_sdm %>%
  filter(!is.na(rcp85)) %>%
  mutate(category = factor(rcp85, levels = 0:7),
         scenario = "RCP 8.5") %>%
  select(category, scenario, CH_no_habitat)

plot_all <- bind_rows(long_all_45, long_all_85)

p_all <- ggplot(plot_all, aes(x = category, y = CH_no_habitat, fill = scenario)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, outlier.size = 0.3, fill = "white", alpha = 0.6) +
  facet_wrap(~ scenario, ncol = 1) +
  scale_fill_manual(values = c("RCP 4.5" = "steelblue", "RCP 8.5" = "firebrick")) +
  labs(title = "All grassland species combined",
       subtitle = paste0(nrow(all_sdm), " route-species observations across ",
                         length(unique(all_sdm$route)), " unique routes and ",
                         length(unique(all_sdm$species_code)), " species"),
       x = "SDM classified change category",
       y = "Residual (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

overall_plot_file <- file.path(plot_dir, "ALL_species_CH_no_habitat_by_rcp.png")
ggsave(overall_plot_file, p_all, width = 8, height = 8, dpi = 150)
cat("Saved overall plot:", overall_plot_file, "\n")

# # Also save the combined CSV
# write.csv(all_sdm,
#           here::here("output", "species_routes_sdm", "ALL_species_route_CH_sdm.csv"),
#           row.names = FALSE)
# cat("Saved combined CSV: output/species_routes_sdm/ALL_species_route_CH_sdm.csv\n")

# ==========================================================================
# Repeat plots using CH (Full with habitat change)
# ==========================================================================

cat("\nGenerating CH violin plots (Full with habitat change)...\n")

# Per-species CH violin plots
for (sf in sdm_files) {
  sp_data <- read.csv(sf)

  abbr <- unique(sp_data$species_code)
  sp_name <- unique(sp_data$species)

  # Skip if plot already exists
  plot_file <- file.path(plot_dir, paste0(abbr, "_CH_by_rcp.png"))
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
    select(category, scenario, CH)

  long_85 <- sp_data %>%
    filter(!is.na(rcp85)) %>%
    mutate(category = factor(rcp85, levels = 0:7),
           scenario = "RCP 8.5") %>%
    select(category, scenario, CH)

  plot_data <- bind_rows(long_45, long_85)

  # Create faceted violin plot
  p <- ggplot(plot_data, aes(x = category, y = CH, fill = scenario)) +
    geom_violin(trim = FALSE, scale = "width") +
    geom_boxplot(width = 0.1, outlier.size = 0.5, fill = "white", alpha = 0.6) +
    facet_wrap(~ scenario, ncol = 1) +
    scale_fill_manual(values = c("RCP 4.5" = "steelblue", "RCP 8.5" = "firebrick")) +
    labs(title = paste0(sp_name, " (", abbr, ")"),
         subtitle = paste0("n = ", nrow(sp_data), " routes"),
         x = "SDM classified change category",
         y = "Full with habitat change (%)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  ggsave(plot_file, p, width = 8, height = 8, dpi = 150)
  cat("  Saved:", basename(plot_file), "\n")
}

# Overall CH violin plot (all species combined)
long_all_45_ch <- all_sdm %>%
  filter(!is.na(rcp45)) %>%
  mutate(category = factor(rcp45, levels = 0:7),
         scenario = "RCP 4.5") %>%
  select(category, scenario, CH)

long_all_85_ch <- all_sdm %>%
  filter(!is.na(rcp85)) %>%
  mutate(category = factor(rcp85, levels = 0:7),
         scenario = "RCP 8.5") %>%
  select(category, scenario, CH)

plot_all_ch <- bind_rows(long_all_45_ch, long_all_85_ch)

p_all_ch <- ggplot(plot_all_ch, aes(x = category, y = CH, fill = scenario)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, outlier.size = 0.3, fill = "white", alpha = 0.6) +
  facet_wrap(~ scenario, ncol = 1) +
  scale_fill_manual(values = c("RCP 4.5" = "steelblue", "RCP 8.5" = "firebrick")) +
  labs(title = "All grassland species combined",
       subtitle = paste0(nrow(all_sdm), " route-species observations across ",
                         length(unique(all_sdm$route)), " unique routes and ",
                         length(unique(all_sdm$species_code)), " species"),
       x = "SDM classified change category",
       y = "Full with habitat change (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

overall_ch_plot_file <- file.path(plot_dir, "ALL_species_CH_by_rcp.png")
ggsave(overall_ch_plot_file, p_all_ch, width = 8, height = 8, dpi = 150)
cat("Saved overall CH plot:", overall_ch_plot_file, "\n")

# Overall CH_dif violin plot (all species combined) -------------------------
long_all_45_dif <- all_sdm %>%
  filter(!is.na(rcp45)) %>%
  mutate(category = factor(rcp45, levels = 0:7),
         scenario = "RCP 4.5") %>%
  select(category, scenario, CH_dif)

long_all_85_dif <- all_sdm %>%
  filter(!is.na(rcp85)) %>%
  mutate(category = factor(rcp85, levels = 0:7),
         scenario = "RCP 8.5") %>%
  select(category, scenario, CH_dif)

plot_all_dif <- bind_rows(long_all_45_dif, long_all_85_dif)

p_all_dif <- ggplot(plot_all_dif, aes(x = category, y = CH_dif, fill = scenario)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, outlier.size = 0.3, fill = "white", alpha = 0.6) +
  facet_wrap(~ scenario, ncol = 1) +
  scale_fill_manual(values = c("RCP 4.5" = "steelblue", "RCP 8.5" = "firebrick")) +
  labs(title = "All grassland species combined",
       subtitle = paste0(nrow(all_sdm), " route-species observations across ",
                         length(unique(all_sdm$route)), " unique routes and ",
                         length(unique(all_sdm$species_code)), " species"),
       x = "SDM classified change category",
       y = "CH difference (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

overall_dif_plot_file <- file.path(plot_dir, "ALL_species_CH_dif_by_rcp.png")
ggsave(overall_dif_plot_file, p_all_dif, width = 8, height = 8, dpi = 150)
cat("Saved overall CH_dif plot:", overall_dif_plot_file, "\n")

# ==========================================================================
# Statistical analysis: differences between SDM change categories (combined)
# ==========================================================================

cat("\n=== Statistical analysis (combined across all species) ===\n")

stats_dir <- here::here("output", "species_routes_sdm_stats")
if (!dir.exists(stats_dir)) dir.create(stats_dir, recursive = TRUE)

stats_file <- file.path(stats_dir, "category_stats.txt")
sink(stats_file)

# Helper: run Kruskal-Wallis + pairwise Wilcoxon (BH-corrected) -------------
run_category_tests <- function(df, response, group, scenario_label, response_label) {
  d <- df %>%
    filter(!is.na(.data[[response]]), !is.na(.data[[group]])) %>%
    mutate(grp = factor(.data[[group]]))

  cat("\n--------------------------------------------------\n")
  cat(scenario_label, "|", response_label, "by", group, "\n")
  cat("--------------------------------------------------\n")

  # Need at least 2 groups with data
  if (length(unique(d$grp)) < 2) {
    cat("  Not enough groups for testing.\n")
    return(invisible(NULL))
  }

  # Group sizes and medians
  summ_tab <- d %>%
    group_by(grp) %>%
    summarise(n = n(),
              median = median(.data[[response]]),
              mean = mean(.data[[response]]),
              .groups = "drop")
  cat("\nGroup summary:\n")
  print(as.data.frame(summ_tab), row.names = FALSE)

  # Kruskal-Wallis omnibus test
  kw <- kruskal.test(d[[response]] ~ d$grp)
  cat("\nKruskal-Wallis test:\n")
  cat("  chi-squared =", round(kw$statistic, 3),
      ", df =", kw$parameter,
      ", p-value =", format.pval(kw$p.value, digits = 4), "\n")

  # Pairwise Wilcoxon (BH-corrected) if omnibus is informative
  pw <- pairwise.wilcox.test(d[[response]], d$grp, p.adjust.method = "BH")
  cat("\nPairwise Wilcoxon (BH-adjusted p-values):\n")
  print(round(pw$p.value, 4))

  invisible(list(kruskal = kw, pairwise = pw))
}

# Build combined long data with category as numeric per scenario -----------
analysis_45 <- all_sdm %>%
  filter(!is.na(rcp45)) %>%
  transmute(category = rcp45, CH, CH_no_habitat, CH_dif)

analysis_85 <- all_sdm %>%
  filter(!is.na(rcp85)) %>%
  transmute(category = rcp85, CH, CH_no_habitat, CH_dif)

# Run tests across all 8 categories ----------------------------------------
cat("\n##################################################\n")
cat("# PART 1: All 8 SDM change categories (0-7)\n")
cat("##################################################\n")

for (resp in c("CH", "CH_no_habitat", "CH_dif")) {
  resp_label <- switch(resp,
                       CH = "Full with habitat change",
                       CH_no_habitat = "Residual",
                       CH_dif = "CH difference")
  run_category_tests(analysis_45, resp, "category", "RCP 4.5", resp_label)
  run_category_tests(analysis_85, resp, "category", "RCP 8.5", resp_label)
}

# ==========================================================================
# PART 2: Grouped categories — contraction / stable / expansion
#   1,2,3 = Contraction ; 4 = Stable ; 5,6,7 = Expansion
#   (category 0 = never suitable, excluded)
# ==========================================================================

cat("\n##################################################\n")
cat("# PART 2: Grouped categories (Contraction/Stable/Expansion)\n")
cat("#   Contraction = 1,2,3 | Stable = 4 | Expansion = 5,6,7\n")
cat("#   (category 0 = never suitable, excluded)\n")
cat("##################################################\n")

group_category <- function(x) {
  dplyr::case_when(
    x %in% c(1, 2, 3) ~ "Contraction",
    x == 4            ~ "Stable",
    x %in% c(5, 6, 7) ~ "Expansion",
    TRUE              ~ NA_character_
  )
}

analysis_45_grp <- analysis_45 %>%
  mutate(change_group = factor(group_category(category),
                               levels = c("Contraction", "Stable", "Expansion")))

analysis_85_grp <- analysis_85 %>%
  mutate(change_group = factor(group_category(category),
                               levels = c("Contraction", "Stable", "Expansion")))

for (resp in c("CH", "CH_no_habitat", "CH_dif")) {
  resp_label <- switch(resp,
                       CH = "Full with habitat change",
                       CH_no_habitat = "Residual",
                       CH_dif = "CH difference")
  run_category_tests(analysis_45_grp, resp, "change_group", "RCP 4.5", resp_label)
  run_category_tests(analysis_85_grp, resp, "change_group", "RCP 8.5", resp_label)
}

sink()
cat("Saved statistics to:", stats_file, "\n")

# ==========================================================================
# Grouped violin plots (Contraction / Stable / Expansion)
# ==========================================================================

make_grouped_violin <- function(df45, df85, response, response_label, file_suffix) {
  long45 <- df45 %>%
    filter(!is.na(change_group)) %>%
    transmute(change_group, scenario = "RCP 4.5", value = .data[[response]])
  long85 <- df85 %>%
    filter(!is.na(change_group)) %>%
    transmute(change_group, scenario = "RCP 8.5", value = .data[[response]])
  pdata <- bind_rows(long45, long85)

  p <- ggplot(pdata, aes(x = change_group, y = value, fill = scenario)) +
    geom_violin(trim = FALSE, scale = "width") +
    geom_boxplot(width = 0.1, outlier.size = 0.3, fill = "white", alpha = 0.6) +
    facet_wrap(~ scenario, ncol = 1) +
    scale_fill_manual(values = c("RCP 4.5" = "steelblue", "RCP 8.5" = "firebrick")) +
    labs(title = "All grassland species combined",
         subtitle = paste0(nrow(all_sdm), " route-species observations across ",
                           length(unique(all_sdm$route)), " unique routes and ",
                           length(unique(all_sdm$species_code)), " species"),
         x = "Range change group",
         y = paste0(response_label, " (%)")) +
    theme_minimal() +
    theme(legend.position = "none")

  out_f <- file.path(plot_dir, paste0("ALL_species_", file_suffix, "_by_group.png"))
  ggsave(out_f, p, width = 7, height = 8, dpi = 150)
  cat("Saved grouped plot:", out_f, "\n")
}

make_grouped_violin(analysis_45_grp, analysis_85_grp, "CH_no_habitat", "Residual", "CH_no_habitat")
make_grouped_violin(analysis_45_grp, analysis_85_grp, "CH", "Full with habitat change", "CH")
make_grouped_violin(analysis_45_grp, analysis_85_grp, "CH_dif", "CH difference", "CH_dif")

cat("\n=== Statistical analysis complete ===\n")

