# 3_statistical_analysis.R

library(here)
library(tidyverse)

here::i_am("code/3_statistical_analysis.R")

out_dir <- here::here("output", "species_routes_sdm")
# Read all SDM CSVs
sdm_files <- list.files(out_dir, pattern = "_route_CH_sdm\\.csv$",
                        full.names = TRUE)

cat("\nCombining all species SDM files...\n")

all_sdm <- sdm_files %>%
  map_dfr(read.csv)

cat("Total rows:", nrow(all_sdm), "\n")
cat("Species:", length(unique(all_sdm$species_code)), "\n")

# ==========================================================================
# Statistical analysis: differences between SDM change categories
# Focus response: CH_no_habitat (Residual); CH (Full) included for reference
# ==========================================================================

stats_dir <- here::here("output", "species_routes_sdm_stats")
if (!dir.exists(stats_dir)) dir.create(stats_dir, recursive = TRUE)

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

  # Group sizes and summary
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
  cat("\nKruskal-Wallis: chi-sq =", round(kw$statistic, 2),
      ", df =", kw$parameter,
      ", p =", format.pval(kw$p.value, digits = 3), "\n")

  # Pairwise Wilcoxon (BH-corrected)
  pw <- pairwise.wilcox.test(d[[response]], d$grp, p.adjust.method = "BH")
  cat("\nPairwise Wilcoxon (BH-adjusted p-values):\n")
  print(round(pw$p.value, 4))

  invisible(list(kruskal = kw, pairwise = pw))
}

# Build combined analysis data ---------------------------------------------
analysis_45 <- all_sdm %>%
  filter(!is.na(rcp45)) %>%
  transmute(route, species_code, category = rcp45, CH, CH_no_habitat, CH_dif,
            route_num = as.integer(sub("^\\d+-", "", route)))

analysis_85 <- all_sdm %>%
  filter(!is.na(rcp85)) %>%
  transmute(route, species_code, category = rcp85, CH, CH_no_habitat, CH_dif,
            route_num = as.integer(sub("^\\d+-", "", route)))

species_list <- unique(all_sdm$species_code)


# --------------------------------------------------------------------------
# PART 1: All 8 SDM categories (0-7) — all species ####
# --------------------------------------------------------------------------

stats_file_8cat <- file.path(stats_dir, "category_stats_8categories.txt")
sink(stats_file_8cat)

cat("\n##################################################\n")
cat("# PART 1: All 8 SDM categories (0-7) — all species\n")
cat("#   Response: CH_no_habitat (Residual) and CH (Full)\n")
cat("##################################################\n")

for (resp in c("CH_no_habitat", "CH")) {
  resp_label <- switch(resp,
                       CH = "Full with habitat change",
                       CH_no_habitat = "Residual")
  run_category_tests(analysis_45, resp, "category", "RCP 4.5", resp_label)
  run_category_tests(analysis_85, resp, "category", "RCP 8.5", resp_label)
}

sink()
cat("Saved Part 1 (all-species 8-category) statistics to:",
    stats_file_8cat, "\n")


# --------------------------------------------------------------------------
# PART 2: All 8 SDM categories (0-7) — per species ####
# --------------------------------------------------------------------------

stats_file_species_8cat <- file.path(stats_dir,
                                      "category_stats_by_species_8categories.txt")
sink(stats_file_species_8cat)

cat("\n##################################################\n")
cat("# PART 2: All 8 SDM categories (0-7) — per species\n")
cat("#   Response: CH_no_habitat (Residual) and CH (Full)\n")
cat("##################################################\n")

for (sp in species_list) {

  sp_name <- unique(all_sdm$species[all_sdm$species_code == sp])

  cat("\n==========================================================\n")
  cat(" ", sp_name, " (", sp, ")\n")
  cat("==========================================================\n")

  # Species-specific data
  sp_45 <- analysis_45 %>% filter(species_code == sp)
  sp_85 <- analysis_85 %>% filter(species_code == sp)

  cat("\n--- RCP 4.5: All 8 categories ---\n")
  run_category_tests(sp_45, "CH_no_habitat", "category",
                     paste0(sp, " | RCP 4.5"), "Residual")
  run_category_tests(sp_45, "CH", "category",
                     paste0(sp, " | RCP 4.5"), "Full with habitat change")

  cat("\n--- RCP 8.5: All 8 categories ---\n")
  run_category_tests(sp_85, "CH_no_habitat", "category",
                     paste0(sp, " | RCP 8.5"), "Residual")
  run_category_tests(sp_85, "CH", "category",
                     paste0(sp, " | RCP 8.5"), "Full with habitat change")
}

sink()
cat("Saved Part 2 (per-species 8-category) statistics to:",
    stats_file_species_8cat, "\n")


# --------------------------------------------------------------------------
# PART 3: Create grouped categories ####
#   Contraction = 1,2,3 | Stable = 4 | Expansion = 5,6,7
#   (category 0 = never suitable, excluded)
# --------------------------------------------------------------------------

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


# --------------------------------------------------------------------------
# PART 4: Grouped categories (Contraction/Stable/Expansion) — all species ####
# --------------------------------------------------------------------------

stats_file_grp <- file.path(stats_dir, "category_stats_grouped.txt")
sink(stats_file_grp)

cat("\n##################################################\n")
cat("# PART 4: Grouped categories — all species\n")
cat("#   Contraction = 1,2,3 | Stable = 4 | Expansion = 5,6,7\n")
cat("#   (category 0 = never suitable, excluded)\n")
cat("#   Response: CH_no_habitat (Residual) and CH (Full)\n")
cat("##################################################\n")

for (resp in c("CH_no_habitat", "CH")) {
  resp_label <- switch(resp,
                       CH = "Full with habitat change",
                       CH_no_habitat = "Residual")
  run_category_tests(analysis_45_grp, resp, "change_group", "RCP 4.5", resp_label)
  run_category_tests(analysis_85_grp, resp, "change_group", "RCP 8.5", resp_label)
}

sink()
cat("Saved Part 4 (all-species grouped) statistics to:",
    stats_file_grp, "\n")


# --------------------------------------------------------------------------
# PART 5: Grouped categories (Contraction/Stable/Expansion) — per species ####
# --------------------------------------------------------------------------

stats_file_species_grp <- file.path(stats_dir,
                                     "category_stats_by_species_grouped.txt")
sink(stats_file_species_grp)

cat("\n##################################################\n")
cat("# PART 5: Grouped categories — per species\n")
cat("#   Contraction = 1,2,3 | Stable = 4 | Expansion = 5,6,7\n")
cat("#   Response: CH_no_habitat (Residual) and CH (Full)\n")
cat("##################################################\n")

for (sp in species_list) {

  sp_name <- unique(all_sdm$species[all_sdm$species_code == sp])

  cat("\n==========================================================\n")
  cat(" ", sp_name, " (", sp, ")\n")
  cat("==========================================================\n")

  sp_45 <- analysis_45_grp %>% filter(species_code == sp)
  sp_85 <- analysis_85_grp %>% filter(species_code == sp)

  cat("\n--- RCP 4.5: Grouped categories ---\n")
  run_category_tests(sp_45, "CH_no_habitat", "change_group",
                     paste0(sp, " | RCP 4.5"), "Residual")
  run_category_tests(sp_45, "CH", "change_group",
                     paste0(sp, " | RCP 4.5"), "Full with habitat change")

  cat("\n--- RCP 8.5: Grouped categories ---\n")
  run_category_tests(sp_85, "CH_no_habitat", "change_group",
                     paste0(sp, " | RCP 8.5"), "Residual")
  run_category_tests(sp_85, "CH", "change_group",
                     paste0(sp, " | RCP 8.5"), "Full with habitat change")
}

sink()
cat("Saved Part 5 (per-species grouped) statistics to:",
    stats_file_species_grp, "\n")

cat("\n=== Statistical analysis complete ===\n")
