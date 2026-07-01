# End-to-end smoke test for the scenario-aware workflow:
#   3_assign_model_scenarios.R -> 4_scenario_fit.R -> 5_scenario_diagnostics.R
#
# This script runs the scenario-aware part of the pipeline on a small,
# diagnostic-rich subset before launching the full species group. It still fits
# real Stan models, so runtime can be substantial for the larger species in the
# subset.
#
# NOTE ON BASELINE DIAGNOSTICS:
# 3_assign_model_scenarios.R chooses each species' model scenario from baseline
# diagnostics produced by 1_baseline_NB_fit.R -> 2_baseline_diagnostics.R (in
# output/model_diagnostics/). If those are absent, 3_assign_model_scenarios.R
# degrades gracefully and assigns the default NB_test scenario to every species,
# so this test still runs end-to-end but only exercises the NB path. To test
# real scenario selection, run the baseline phase first:
#   source("code/1_baseline_NB_fit.R")        # baseline NB fit (all species)
#   source("code/2_baseline_diagnostics.R")   # writes output/model_diagnostics/

library(here)

here::i_am("code/tests/test_scenario_pipeline_subset.R")

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

cat("\n=== Step 3: prepare model-scenario assignments ===\n")
source(here::here("code", "3_assign_model_scenarios.R"), local = FALSE)

cat("\n=== Step 4: fit selected species and write route summaries ===\n")
cat("Subset species codes:", paste(only_species_codes, collapse = ", "), "\n")
source(here::here("code", "4_scenario_fit.R"), local = FALSE)

cat("\n=== Step 5: diagnose selected scenario fits ===\n")
source(here::here("code", "5_scenario_diagnostics.R"), local = FALSE)

cat("\n=== End-to-end subset test complete ===\n")
