# End-to-end smoke test for the scenario-aware workflow:
#   1a_prepare_model_scenarios.R -> 1_CH_no_habitat_routes.R ->
#   1b_model_diagnostics_by_scenario.R
#
# This script runs the full pipeline on a small, diagnostic-rich subset before
# launching the full species group. It still fits real Stan models, so runtime
# can be substantial for the larger species in the subset.

library(here)

here::i_am("code/test_1a_1_1b_subset.R")

# Around 10 species spanning clean, marginal, excess-zero, and difficult cases.
# Edit this vector if you want a faster or harder test.
only_species_codes <- c(
  "BOSP", # clean control
  "LEPC", # clean low/moderate count
  "WTHA", # clean rare species
  "CCLO", # route-ZINB candidate
  "UPSA", # difficult route+abundance ZINB candidate
  "BOBO", # difficult route+abundance ZINB candidate
  "CASP", # excess-zero candidate
  "DICK", # excess-zero candidate
  "LARB", # poor convergence / excess-zero candidate
  "SWHA"  # opposite-direction PPC spread issue candidate
)

cat("\n=== Step 1a: prepare model-scenario assignments ===\n")
source(here::here("code", "1a_prepare_model_scenarios.R"), local = FALSE)

cat("\n=== Step 1: fit selected species and write route summaries ===\n")
cat("Subset species codes:", paste(only_species_codes, collapse = ", "), "\n")
source(here::here("code", "1_CH_no_habitat_routes.R"), local = FALSE)

cat("\n=== Step 1b: diagnose selected scenario fits ===\n")
source(here::here("code", "1b_model_diagnostics_by_scenario.R"), local = FALSE)

cat("\n=== End-to-end subset test complete ===\n")
