# Compare population-change trend summaries across the three fitted models.
# Reads the three *_summ_fit.rds files in output/ and produces a comparison
# table of the key trend quantities (T, T_no_habitat, T_dif, CH, CH_no_habitat,
# CH_dif) on the same scale.

library(tidyverse)
library(here)
library(knitr)

here::i_am("code/2_compare_trends.R")

# List and name the three fits
fits <- c(
  "Developed"       = "Bairds_Sparrow_developed_2010_2024_summ_fit.rds",
  "Grassland"       = "Bairds_Sparrow_grassland_2010_2024_summ_fit.rds",
  "iCAR (no cov.)"  = "Bairds_Sparrow_iCAR_New_2010_2024_summ_fit.rds"
)

summ_list <- lapply(fits, function(f) readRDS(here::here("output", f)))

# Show the parameter names each file actually has for the trend quantities, so
# we can build the table robustly even if the iCAR-only model has fewer vars.
target_vars <- c("T", "T_no_habitat", "T_dif",
                 "CH", "CH_no_habitat", "CH_dif",
                 "BETA", "ALPHA",
                 "rho_BETA_hab", "rho_ALPHA_hab")

fmt <- function(x, digits = 3) formatC(x, format = "f", digits = digits)

extract_row <- function(s, model_name) {
  s %>%
    filter(variable %in% target_vars) %>%
    mutate(ci90 = paste0("[", fmt(q5), ", ", fmt(q95), "]"),
           model = model_name) %>%
    select(model, variable, mean, median, sd, q5, q95, rhat, ess_bulk)
}

comparison <- purrr::imap_dfr(summ_list, extract_row)

# Wide view: rows = variable, cols = model, cell = "mean [q5, q95]"
wide <- comparison %>%
  mutate(cell = paste0(fmt(mean), " [", fmt(q5), ", ", fmt(q95), "]")) %>%
  select(variable, model, cell) %>%
  pivot_wider(names_from = model, values_from = cell) %>%
  # keep variables in a consistent, readable order
  mutate(variable = factor(variable, levels = target_vars)) %>%
  arrange(variable)

cat("\n=== Population trend comparison (mean [90% CI]) ===\n")
print(kable(wide, format = "pipe"))

# Long view with diagnostics, for any model that reports the row
cat("\n=== Long form with rhat / ESS ===\n")
print(kable(comparison %>%
              mutate(across(c(mean, median, sd, q5, q95), ~fmt(.x))) %>%
              mutate(rhat = fmt(rhat, 3),
                     ess_bulk = round(ess_bulk)),
            format = "pipe"))

# Also write a markdown file for the record
out_md <- here::here("output", "trend_comparison.md")
lines <- c(
  "# Population trend comparison — Baird's Sparrow (2010–2024)",
  "",
  "Three route-level iCAR slope-model fits, compared on the trend quantities",
  "reported by each fit.",
  "",
  "- **Developed**: single-covariate model, proportion developed land.",
  "- **Grassland**: single-covariate model, proportion grassland.",
  "- **iCAR (no cov.)**: null spatial model, no habitat covariates.",
  "",
  "## Trend quantities (mean [90% CI])",
  "",
  paste(kable(wide, format = "pipe"), collapse = "\n"),
  "",
  "## All rows with diagnostics",
  "",
  paste(kable(comparison %>%
                mutate(across(c(mean, median, sd, q5, q95), ~fmt(.x))) %>%
                mutate(rhat = fmt(rhat, 3),
                       ess_bulk = round(ess_bulk)),
              format = "pipe"),
        collapse = "\n"),
  ""
)
writeLines(lines, out_md)
cat("\nWrote ", out_md, "\n", sep = "")
