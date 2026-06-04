# Per-route trend comparison across the three fitted models.
# Builds two tables (T and CH) with routes as rows and the following columns:
#   T_null, T_dev, T_no_hab_dev, T_dif_dev, T_gr, T_no_hab_gr, T_dif_gr
#   CH_null, CH_dev, CH_no_hab_dev, CH_dif_dev, CH_gr, CH_no_hab_gr, CH_dif_gr
#
# Route-level trend is computed from the posterior mean of beta[s] / beta_resid[s]:
#   T_route[s]  = 100 * (exp(beta[s])       - 1)   (%/yr)
#   CH_route[s] = 100 * (exp(beta[s] * dt)  - 1)   (% over the time window)
# where dt = nyears - 1 (the length of the window in years).
#
# T_no_habitat_route uses beta_resid[s] instead of beta[s]; T_dif is the diff.

library(tidyverse)
library(here)
library(knitr)

here::i_am("code/3_compare_route_trends.R")

fits <- list(
  null = list(
    label = "null",
    file  = "Bairds_Sparrow_iCAR_New_2010_2024_summ_fit.rds"
  ),
  dev = list(
    label = "dev",
    file  = "Bairds_Sparrow_developed_2010_2024_summ_fit.rds"
  ),
  gr = list(
    label = "gr",
    file  = "Bairds_Sparrow_grassland_2010_2024_summ_fit.rds"
  )
)

# length of the time window (in years). 2010..2024 inclusive = 15 years, so dt = 14.
firstYear <- 2010
lastYear  <- 2024
dt <- lastYear - firstYear   # 14

# helpers
route_idx <- function(var) as.integer(str_extract(var, "\\d+"))

# Extract route-level beta / beta_resid posterior means for one fit.
# Returns a long tibble with columns: model, routeF, beta, beta_resid (may be NA).
extract_route_betas <- function(summ, model_label) {
  b <- summ %>%
    filter(str_detect(variable, "^beta\\[")) %>%
    transmute(routeF = route_idx(variable), beta = mean)

  br <- summ %>%
    filter(str_detect(variable, "^beta_resid\\[")) %>%
    transmute(routeF = route_idx(variable), beta_resid = mean)

  df <- full_join(b, br, by = "routeF") %>%
    mutate(model = model_label)

  df
}

per_route <- purrr::imap_dfr(fits, function(fs, key) {
  path <- here::here("output", fs$file)
  if (!file.exists(path)) {
    warning("Missing fit file: ", path)
    return(tibble())
  }
  summ <- readRDS(path)
  extract_route_betas(summ, fs$label)
})

# Convert beta / beta_resid to T and CH quantities
pct_yr  <- function(b) 100 * (exp(b) - 1)
pct_cum <- function(b) 100 * (exp(b * dt) - 1)

per_route <- per_route %>%
  mutate(T_route          = pct_yr(beta),
         T_no_hab_route   = pct_yr(beta_resid),
         T_dif_route      = T_route - T_no_hab_route,
         CH_route         = pct_cum(beta),
         CH_no_hab_route  = pct_cum(beta_resid),
         CH_dif_route     = CH_route - CH_no_hab_route)

# Pivot to wide: one row per routeF, columns per model
T_table <- per_route %>%
  select(routeF, model, T_route, T_no_hab_route, T_dif_route) %>%
  pivot_wider(id_cols = routeF, names_from = model,
              values_from = c(T_route, T_no_hab_route, T_dif_route),
              names_glue = "{.value}_{model}") %>%
  arrange(routeF) %>%
  # final column order: null first, then dev triad, then gr triad
  select(routeF,
         T_null        = T_route_null,
         T_dev         = T_route_dev,
         T_no_hab_dev  = T_no_hab_route_dev,
         T_dif_dev     = T_dif_route_dev,
         T_gr          = T_route_gr,
         T_no_hab_gr   = T_no_hab_route_gr,
         T_dif_gr      = T_dif_route_gr)

CH_table <- per_route %>%
  select(routeF, model, CH_route, CH_no_hab_route, CH_dif_route) %>%
  pivot_wider(id_cols = routeF, names_from = model,
              values_from = c(CH_route, CH_no_hab_route, CH_dif_route),
              names_glue = "{.value}_{model}") %>%
  arrange(routeF) %>%
  select(routeF,
         CH_null        = CH_route_null,
         CH_dev         = CH_route_dev,
         CH_no_hab_dev  = CH_no_hab_route_dev,
         CH_dif_dev     = CH_dif_route_dev,
         CH_gr          = CH_route_gr,
         CH_no_hab_gr   = CH_no_hab_route_gr,
         CH_dif_gr      = CH_dif_route_gr)

fmt <- function(x, digits = 3) ifelse(is.na(x), "–", formatC(x, format = "f", digits = digits))

T_fmt  <- T_table  %>% mutate(across(-routeF, fmt))
CH_fmt <- CH_table %>% mutate(across(-routeF, fmt))

cat("\n=== Per-route T (%/yr) ===\n")
print(kable(T_fmt, format = "pipe"))

cat("\n=== Per-route CH (% over ", dt, " yr) ===\n", sep = "")
print(kable(CH_fmt, format = "pipe"))

# Write a markdown file
out_md <- here::here("output", "route_trend_comparison.md")
lines <- c(
  "# Per-route trend comparison — Baird's Sparrow (2010–2024)",
  "",
  "Route-level trends derived from the posterior mean of `beta[s]` and",
  "`beta_resid[s]` for each fit:",
  "",
  "- `T_route[s]   = 100 * (exp(beta[s]) - 1)`   (%/yr)",
  paste0("- `CH_route[s]  = 100 * (exp(beta[s] * ", dt, ") - 1)`  (% cumulative over the window)"),
  "- `T_no_hab_route[s]` uses `beta_resid[s]` instead of `beta[s]`.",
  "- `T_dif_route[s] = T_route[s] - T_no_hab_route[s]`.",
  "",
  "\"–\" means the quantity is not defined for that model (e.g. the null iCAR",
  "model has no habitat decomposition, so no `no_hab` / `dif` columns).",
  "",
  "## T — annual trend (%/yr)",
  "",
  paste(kable(T_fmt, format = "pipe"), collapse = "\n"),
  "",
  paste0("## CH — cumulative change (% over ", dt, " yr)"),
  "",
  paste(kable(CH_fmt, format = "pipe"), collapse = "\n"),
  ""
)
writeLines(lines, out_md)
cat("\nWrote ", out_md, "\n", sep = "")

# Also save the unformatted tables as CSV for downstream use
write.csv(T_table,  here::here("output", "route_T_comparison.csv"),  row.names = FALSE)
write.csv(CH_table, here::here("output", "route_CH_comparison.csv"), row.names = FALSE)
