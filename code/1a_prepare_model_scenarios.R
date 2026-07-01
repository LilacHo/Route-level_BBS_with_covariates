# Prepare species-to-model scenario assignments for
# code/1_CH_no_habitat_routes.R.
#
# This script turns baseline diagnostics into a reproducible scenario table.
# It is intentionally rule-based and transparent: each species receives a
# model_scenario plus an assignment_reason that should be reviewed before large
# production runs.

library(here)
library(tidyverse)

here::i_am("code/1a_prepare_model_scenarios.R")

# Settings -----------------------------------------------------------------
land_cover <- "grasslands"

rhat_thresh <- 1.01
poor_rhat_thresh <- 1.10
ess_thresh <- 400
low_ess_thresh <- 100
ppc_tail <- 0.05
strong_zero_tail <- 0.01
abundant_mean_thresh <- 5

diag_dir <- here::here("output", "model_diagnostics")
out_dir <- here::here("output", "model_scenario_assignments")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Existing diagnostics were historically written with the singular prefix
# "grassland"; the current land_cover is "grasslands". Try both.
diagnostic_prefixes <- unique(c(land_cover, sub("s$", "", land_cover)))

first_existing <- function(pattern) {
  candidates <- file.path(diag_dir, paste0(diagnostic_prefixes, pattern))
  candidates[file.exists(candidates)][1]
}

conv_file <- first_existing("_convergence_summary.csv")
ppc_file <- first_existing("_posterior_predictive_checks.csv")
sampler_file <- first_existing("_sampler_diagnostics.csv")

if (is.na(conv_file) || is.na(ppc_file)) {
  stop("Could not find baseline convergence and PPC diagnostics in ", diag_dir,
       ". Run diagnostics on the baseline model first.")
}

spp_df <- read.csv(here::here("data", "spp_names_codes_group_aou.csv")) %>%
  filter(Group == land_cover) %>%
  distinct(Common.Name, Code, bbs_english, .keep_all = TRUE)

conv <- read.csv(conv_file) %>%
  transmute(
    species_code,
    baseline_max_rhat = max_rhat,
    baseline_min_ess_bulk = min_ess_bulk,
    baseline_min_ess_tail = min_ess_tail,
    baseline_worst_rhat_param = worst_rhat_param
  )

ppc <- read.csv(ppc_file) %>%
  transmute(
    species_code,
    obs_mean,
    obs_sd,
    obs_prop_zero,
    obs_max,
    p_mean,
    p_sd,
    p_prop_zero,
    p_max
  )

sampler <- if (!is.na(sampler_file)) {
  read.csv(sampler_file) %>%
    transmute(
      species_code,
      baseline_num_divergent = num_divergent,
      baseline_num_max_treedepth = num_max_treedepth,
      baseline_min_ebfmi = min_ebfmi
    )
} else {
  tibble(
    species_code = character(),
    baseline_num_divergent = numeric(),
    baseline_num_max_treedepth = numeric(),
    baseline_min_ebfmi = numeric()
  )
}

classify_scenario <- function(p_prop_zero, p_sd, p_max, p_mean,
                              max_rhat, min_ess_bulk, obs_mean) {
  excess_zero <- !is.na(p_prop_zero) && p_prop_zero < ppc_tail
  strong_excess_zero <- !is.na(p_prop_zero) && p_prop_zero < strong_zero_tail
  abundant <- !is.na(obs_mean) && obs_mean >= abundant_mean_thresh
  high_spread <- (!is.na(p_sd) && p_sd > 1 - ppc_tail) ||
    (!is.na(p_max) && p_max > 1 - ppc_tail)
  low_spread <- (!is.na(p_sd) && p_sd < ppc_tail) ||
    (!is.na(p_max) && p_max < ppc_tail)
  severe_mixing <- (!is.na(max_rhat) && max_rhat > poor_rhat_thresh) ||
    (!is.na(min_ess_bulk) && min_ess_bulk < low_ess_thresh)

  if (strong_excess_zero && high_spread) return("ZINB_route_mu_test")
  if (excess_zero && high_spread && abundant) return("ZINB_route_mu_test")
  if (strong_excess_zero && abundant) return("ZINB_route_test")
  if (low_spread) return("NB_test")
  if (severe_mixing) return("NB_test")
  "NB_test"
}

assignment_reason <- function(model_scenario, p_prop_zero, p_sd, p_max,
                              max_rhat, min_ess_bulk) {
  excess_zero <- !is.na(p_prop_zero) && p_prop_zero < ppc_tail
  high_spread <- (!is.na(p_sd) && p_sd > 1 - ppc_tail) ||
    (!is.na(p_max) && p_max > 1 - ppc_tail)
  low_spread <- (!is.na(p_sd) && p_sd < ppc_tail) ||
    (!is.na(p_max) && p_max < ppc_tail)
  severe_mixing <- (!is.na(max_rhat) && max_rhat > poor_rhat_thresh) ||
    (!is.na(min_ess_bulk) && min_ess_bulk < low_ess_thresh)

  case_when(
    model_scenario == "ZINB_route_mu_test" ~
      "Baseline PPC shows strong excess zeros plus inflated spread/max; use route + abundance zero inflation.",
    model_scenario == "ZINB_route_test" ~
      "Baseline PPC shows strong excess zeros in an abundant species without strong spread/max inflation; use route-level zero inflation first.",
    low_spread ~
      "Baseline PPC suggests under-dispersion or extreme-route issue rather than excess zeros; keep NB and inspect separately.",
    severe_mixing ~
      "Baseline mixing was poor, but PPC did not indicate excess-zero structure; use improved NB settings first.",
    TRUE ~
      "Baseline diagnostics do not require zero inflation; use improved NB model."
  )
}

assignments <- spp_df %>%
  left_join(conv, by = c("Code" = "species_code")) %>%
  left_join(ppc, by = c("Code" = "species_code")) %>%
  left_join(sampler, by = c("Code" = "species_code")) %>%
  rowwise() %>%
  mutate(
    model_scenario = classify_scenario(
      p_prop_zero, p_sd, p_max, p_mean,
      baseline_max_rhat, baseline_min_ess_bulk, obs_mean
    ),
    assignment_reason = assignment_reason(
      model_scenario, p_prop_zero, p_sd, p_max,
      baseline_max_rhat, baseline_min_ess_bulk
    ),
    diagnostic_flags = paste(
      c(
        if (!is.na(baseline_max_rhat) && baseline_max_rhat > rhat_thresh) "rhat",
        if (!is.na(baseline_min_ess_bulk) && baseline_min_ess_bulk < ess_thresh) "ess",
        if (!is.na(p_prop_zero) && p_prop_zero < ppc_tail) "excess_zeros",
        if (!is.na(p_sd) && (p_sd < ppc_tail || p_sd > 1 - ppc_tail)) "sd_ppc",
        if (!is.na(p_max) && (p_max < ppc_tail || p_max > 1 - ppc_tail)) "max_ppc"
      ),
      collapse = ";"
    ),
    diagnostic_flags = ifelse(diagnostic_flags == "", "none", diagnostic_flags)
  ) %>%
  ungroup() %>%
  select(
    Common.Name, Code, bbs_english, Group,
    model_scenario, assignment_reason, diagnostic_flags,
    baseline_max_rhat, baseline_min_ess_bulk, baseline_min_ess_tail,
    baseline_worst_rhat_param,
    obs_mean, obs_prop_zero, p_mean, p_sd, p_prop_zero, p_max,
    any_of(c("baseline_num_divergent", "baseline_num_max_treedepth",
             "baseline_min_ebfmi"))
  )

out_file <- file.path(out_dir, paste0(land_cover, "_model_scenario_assignments.csv"))
write.csv(assignments, out_file, row.names = FALSE)

summary_file <- file.path(out_dir, paste0(land_cover, "_model_scenario_summary.txt"))
summary_lines <- c(
  paste("Model scenario assignment summary for", land_cover),
  paste("Convergence diagnostics:", conv_file),
  paste("PPC diagnostics:", ppc_file),
  "",
  "Assignment counts:",
  capture.output(print(assignments %>% count(model_scenario), row.names = FALSE)),
  "",
  "Rules:",
  "- excess zeros + high replicated SD/max -> ZINB_route_mu_test",
  paste0("- strong excess zeros (< ", strong_zero_tail,
         ") in abundant species (mean >= ", abundant_mean_thresh,
         ") -> ZINB_route_test"),
  "- weak/moderate zero flags in low-count species -> NB_test unless spread/max also fails",
  "- low replicated SD/max, or no clear zero problem -> NB_test",
  "",
  "Review species with diagnostic_flags before full production fitting."
)
writeLines(summary_lines, summary_file)

cat("Wrote scenario assignments to:", out_file, "\n")
cat("Wrote summary to:", summary_file, "\n")
print(assignments %>% count(model_scenario))
