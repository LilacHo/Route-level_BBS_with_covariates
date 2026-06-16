# 4_statistical_analysis.R

library(here)
library(tidyverse)

here::i_am("code/4_statistical_analysis.R")

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

# Grouped categories function
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
# PART 1: Grouped categories (Contraction/Stable/Expansion) — all species
# --------------------------------------------------------------------------
cat("\n##################################################\n")
cat("# PART 1: Grouped categories (Contraction/Stable/Expansion)\n")
cat("#   Contraction = 1,2,3 | Stable = 4 | Expansion = 5,6,7\n")
cat("#   (category 0 = never suitable, excluded)\n")
cat("##################################################\n")

for (resp in c("CH_no_habitat", "CH")) {
  resp_label <- switch(resp,
                       CH = "Full with habitat change",
                       CH_no_habitat = "Residual")
  run_category_tests(analysis_45_grp, resp, "change_group", "RCP 4.5", resp_label)
  run_category_tests(analysis_85_grp, resp, "change_group", "RCP 8.5", resp_label)
}

# --------------------------------------------------------------------------
# PART 2: Grouped categories by route < 800 vs >= 800
# --------------------------------------------------------------------------
cat("\n##################################################\n")
cat("# PART 2: Grouped categories by route number\n")
cat("#   Routes < 800 (established) vs Routes >= 800 (newer/random)\n")
cat("##################################################\n")

for (route_subset in c("lt800", "ge800")) {
  if (route_subset == "lt800") {
    sub_45 <- analysis_45_grp %>% filter(route_num < 800)
    sub_85 <- analysis_85_grp %>% filter(route_num < 800)
    subset_label <- "Route < 800"
  } else {
    sub_45 <- analysis_45_grp %>% filter(route_num >= 800)
    sub_85 <- analysis_85_grp %>% filter(route_num >= 800)
    subset_label <- "Route >= 800"
  }
  
  cat("\n==========================================================\n")
  cat("  ", subset_label, " (RCP4.5 n=", nrow(sub_45),
      ", RCP8.5 n=", nrow(sub_85), ")\n")
  cat("==========================================================\n")
  
  for (resp in c("CH_no_habitat", "CH")) {
    resp_label <- switch(resp,
                         CH = "Full with habitat change",
                         CH_no_habitat = "Residual")
    run_category_tests(sub_45, resp, "change_group",
                       paste0("RCP 4.5 | ", subset_label), resp_label)
    run_category_tests(sub_85, resp, "change_group",
                       paste0("RCP 8.5 | ", subset_label), resp_label)
  }
}

sink()
cat("Saved statistics to:", stats_file, "\n")

# ==========================================================================
# PART 3: Per-species statistical analysis (grouped categories)
# ==========================================================================

cat("\nRunning per-species statistical analysis...\n")

stats_file_species <- file.path(stats_dir, "category_stats_by_species.txt")
sink(stats_file_species)

cat("##################################################\n")
cat("# Per-species: Grouped categories (Contraction/Stable/Expansion)\n")
cat("#              + All 8 SDM categories (0-7)\n")
cat("# Response: CH_no_habitat (Residual) and CH (Full)\n")
cat("##################################################\n")

species_list <- unique(all_sdm$species_code)

for (sp in species_list) {
  sp_name <- unique(all_sdm$species[all_sdm$species_code == sp])
  
  cat("\n==========================================================\n")
  cat(" ", sp_name, "(", sp, ")\n")
  cat("==========================================================\n")
  
  sp_45 <- analysis_45_grp %>% filter(species_code == sp)
  sp_85 <- analysis_85_grp %>% filter(species_code == sp)
  
  # All 8 categories
  cat("\n--- All 8 categories ---\n")
  run_category_tests(sp_45, "CH_no_habitat", "category",
                     paste0(sp, " | RCP 4.5"), "Residual")
  run_category_tests(sp_85, "CH_no_habitat", "category",
                     paste0(sp, " | RCP 8.5"), "Residual")
  run_category_tests(sp_45, "CH", "category",
                     paste0(sp, " | RCP 4.5"), "Full with habitat change")
  run_category_tests(sp_85, "CH", "category",
                     paste0(sp, " | RCP 8.5"), "Full with habitat change")
  
  # Grouped categories
  cat("\n--- Grouped (Contraction/Stable/Expansion) ---\n")
  run_category_tests(sp_45, "CH_no_habitat", "change_group",
                     paste0(sp, " | RCP 4.5"), "Residual")
  run_category_tests(sp_85, "CH_no_habitat", "change_group",
                     paste0(sp, " | RCP 8.5"), "Residual")
  run_category_tests(sp_45, "CH", "change_group",
                     paste0(sp, " | RCP 4.5"), "Full with habitat change")
  run_category_tests(sp_85, "CH", "change_group",
                     paste0(sp, " | RCP 8.5"), "Full with habitat change")
}

sink()
cat("Saved per-species statistics to:", stats_file_species, "\n")

# Grouped violin plots broken down by route number -------------------------
make_grouped_violin_route <- function(df45, df85, response, response_label,
                                      file_suffix, route_filter, route_label) {
  long45 <- df45 %>%
    filter(!is.na(change_group)) %>%
    transmute(change_group, scenario = "RCP 4.5", value = .data[[response]])
  long85 <- df85 %>%
    filter(!is.na(change_group)) %>%
    transmute(change_group, scenario = "RCP 8.5", value = .data[[response]])
  pdata <- bind_rows(long45, long85)
  
  if (nrow(pdata) == 0) return(invisible(NULL))
  
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
    labs(title = paste0("All grassland species — ", route_label),
         subtitle = paste0(nrow(pdata) / 2, " route-species observations"),
         x = "Range change group",
         y = paste0(response_label, " (%)")) +
    theme_minimal() +
    theme(legend.position = "none")
  
  out_f <- file.path(plot_dir,
                     paste0("ALL_species_", file_suffix, "_by_group_", route_filter, ".png"))
  ggsave(out_f, p, width = 7, height = 8, dpi = 150)
  cat("Saved grouped plot:", basename(out_f), "\n")
}

# Route < 800
sub_45_lt <- analysis_45_grp %>% filter(route_num < 800)
sub_85_lt <- analysis_85_grp %>% filter(route_num < 800)
make_grouped_violin_route(sub_45_lt, sub_85_lt, "CH_no_habitat", "Residual", "CH_no_habitat", "lt800", "Route < 800")
make_grouped_violin_route(sub_45_lt, sub_85_lt, "CH", "Full with habitat change", "CH", "lt800", "Route < 800")
make_grouped_violin_route(sub_45_lt, sub_85_lt, "CH_dif", "CH difference", "CH_dif", "lt800", "Route < 800")

# Route >= 800
sub_45_ge <- analysis_45_grp %>% filter(route_num >= 800)
sub_85_ge <- analysis_85_grp %>% filter(route_num >= 800)
make_grouped_violin_route(sub_45_ge, sub_85_ge, "CH_no_habitat", "Residual", "CH_no_habitat", "ge800", "Route >= 800")
make_grouped_violin_route(sub_45_ge, sub_85_ge, "CH", "Full with habitat change", "CH", "ge800", "Route >= 800")
make_grouped_violin_route(sub_45_ge, sub_85_ge, "CH_dif", "CH difference", "CH_dif", "ge800", "Route >= 800")


cat("\n=== Statistical analysis complete ===\n")


