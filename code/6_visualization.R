# 6_visualization.R

library(here)
library(tidyverse)

here::i_am("code/6_visualization.R")

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

# Create violin plot output directory
plot_dir <- here::here("output", "species_routes_sdm_plot")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Read all SDM CSVs and generate bar plots
sdm_files <- list.files(out_dir, pattern = "_route_CH_sdm\\.csv$",
                        full.names = TRUE)

cat("\nGenerating violin plots for", length(sdm_files), "species...\n")

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
  count_data <- plot_data %>%
    group_by(category, scenario) %>%
    summarise(n = n(), .groups = "drop")
  
  p <- ggplot(plot_data, aes(x = category, y = CH_no_habitat, fill = scenario)) +
    geom_violin(trim = FALSE, scale = "width") +
    geom_boxplot(width = 0.1, outlier.size = 0.5, fill = "white", alpha = 0.6) +
    geom_text(data = count_data, aes(x = category, y = Inf, label = n),
              vjust = 1.5, size = 2.5, inherit.aes = FALSE) +
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
  count_data <- plot_data %>%
    group_by(category, scenario) %>%
    summarise(n = n(), .groups = "drop")
  
  p <- ggplot(plot_data, aes(x = category, y = CH, fill = scenario)) +
    geom_violin(trim = FALSE, scale = "width") +
    geom_boxplot(width = 0.1, outlier.size = 0.5, fill = "white", alpha = 0.6) +
    geom_text(data = count_data, aes(x = category, y = Inf, label = n),
              vjust = 1.5, size = 2.5, inherit.aes = FALSE) +
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
  geom_text(data = plot_all %>% group_by(category, scenario) %>% summarise(n = n(), .groups = "drop"),
            aes(x = category, y = Inf, label = n),
            vjust = 1.5, size = 2.5, inherit.aes = FALSE) +
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
  geom_text(data = plot_all_ch %>% group_by(category, scenario) %>% summarise(n = n(), .groups = "drop"),
            aes(x = category, y = Inf, label = n),
            vjust = 1.5, size = 2.5, inherit.aes = FALSE) +
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
  geom_text(data = plot_all_dif %>% group_by(category, scenario) %>% summarise(n = n(), .groups = "drop"),
            aes(x = category, y = Inf, label = n),
            vjust = 1.5, size = 2.5, inherit.aes = FALSE) +
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
# Grouped violin plots (Contraction / Stable / Expansion) — all routes
# ==========================================================================

make_grouped_violin <- function(df45, df85, response, response_label, file_suffix) {
  long45 <- df45 %>%
    filter(!is.na(change_group)) %>%
    transmute(change_group, scenario = "RCP 4.5", value = .data[[response]])
  long85 <- df85 %>%
    filter(!is.na(change_group)) %>%
    transmute(change_group, scenario = "RCP 8.5", value = .data[[response]])
  pdata <- bind_rows(long45, long85)
  
  count_data <- pdata %>%
    group_by(change_group, scenario) %>%
    summarise(n = n(), .groups = "drop")
  
  p <- ggplot(pdata, aes(x = change_group, y = value, fill = scenario)) +
    geom_violin(trim = FALSE, scale = "width") +
    geom_boxplot(width = 0.1, outlier.size = 0.3, fill = "white", alpha = 0.6) +
    geom_text(data = count_data, aes(x = change_group, y = Inf, label = n),
              vjust = 1.5, size = 2.5, inherit.aes = FALSE) +
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



