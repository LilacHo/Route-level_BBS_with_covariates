# Diagnostics for scenario-aware route-level BBS fits from
# code/4_scenario_fit.R.
#
# This script reads the scenario manifest created by step 4, then evaluates the
# saved fit for each species/model scenario:
#   1. convergence: R-hat and ESS
#   2. sampler health: divergences, treedepth, E-BFMI
#   3. posterior predictive checks: mean, SD, zero proportion, max
#   4. optional LOO, when the loo package is installed and log_lik is present

library(here)
library(tidyverse)
library(posterior)
library(cmdstanr)

here::i_am("code/5_scenario_diagnostics.R")

# Settings -----------------------------------------------------------------
land_cover <- "grasslands"
firstYear <- 2010
lastYear <- 2024

rhat_thresh <- 1.01
ess_thresh <- 400
ebfmi_thresh <- 0.2
n_ppc_draws <- 500

fit_dir <- here::here("output", "model_fits_by_scenario")
stan_data_dir <- here::here("data", "stan_data_by_scenario")
diag_dir <- here::here("output", "model_diagnostics_by_scenario")
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

manifest_file <- file.path(fit_dir, paste0(land_cover, "_scenario_manifest.csv"))
if (!file.exists(manifest_file)) {
  stop("Scenario manifest not found. Run code/4_scenario_fit.R first: ",
       manifest_file)
}

species_to_f <- function(sp) {
  gsub("'", "", gsub(" ", "_", sp, fixed = TRUE), fixed = TRUE)
}

fit_output_files_exist <- function(fit) {
  files <- tryCatch(
    fit$output_files(include_failed = FALSE),
    error = function(e) character()
  )
  length(files) > 0 && all(file.exists(files))
}

species_base <- function(sp_f, model_scenario, min_max_route_years) {
  paste0(sp_f, "_", land_cover, "_", model_scenario, "_",
         firstYear, "_", lastYear, "_minmax", min_max_route_years)
}

convergence_summary <- function(summ, sp, sp_code, model_scenario) {
  needed <- c("rhat", "ess_bulk", "ess_tail")
  if (!all(needed %in% names(summ))) return(NULL)

  s <- summ %>%
    filter(!is.na(rhat)) %>%
    filter(!str_detect(variable, "^(E|log_lik|zi_obs)\\["))

  tibble(
    species = sp,
    species_code = sp_code,
    model_scenario = model_scenario,
    n_parameters_checked = nrow(s),
    max_rhat = max(s$rhat, na.rm = TRUE),
    n_rhat_gt_thresh = sum(s$rhat > rhat_thresh, na.rm = TRUE),
    prop_rhat_ok = mean(s$rhat <= rhat_thresh, na.rm = TRUE),
    min_ess_bulk = min(s$ess_bulk, na.rm = TRUE),
    min_ess_tail = min(s$ess_tail, na.rm = TRUE),
    n_ess_lt_thresh = sum(pmin(s$ess_bulk, s$ess_tail, na.rm = TRUE) < ess_thresh,
                          na.rm = TRUE),
    worst_rhat_param = s$variable[which.max(s$rhat)]
  )
}

sampler_diagnostics <- function(fit, sp, sp_code, model_scenario) {
  ds <- tryCatch(fit$diagnostic_summary(quiet = TRUE),
                 error = function(e) NULL)
  if (is.null(ds)) return(NULL)

  n_div <- sum(ds$num_divergent)
  n_draws <- tryCatch(posterior::ndraws(fit$draws()), error = function(e) NA_integer_)

  tibble(
    species = sp,
    species_code = sp_code,
    model_scenario = model_scenario,
    num_divergent = n_div,
    divergence_rate = if (!is.na(n_draws) && n_draws > 0) n_div / n_draws else NA_real_,
    num_max_treedepth = sum(ds$num_max_treedepth),
    min_ebfmi = min(ds$ebfmi, na.rm = TRUE),
    n_chains_low_ebfmi = sum(ds$ebfmi < ebfmi_thresh, na.rm = TRUE)
  )
}

posterior_predictive_check <- function(fit, stan_data, sp, sp_code,
                                       model_scenario, sp_f) {
  draws_E <- tryCatch(posterior::as_draws_matrix(fit$draws("E")),
                      error = function(e) NULL)
  draws_phi <- tryCatch(posterior::as_draws_matrix(fit$draws("phi")),
                        error = function(e) NULL)
  if (is.null(draws_E) || is.null(draws_phi)) return(NULL)

  draws_zi <- tryCatch(posterior::as_draws_matrix(fit$draws("zi")),
                       error = function(e) NULL)
  draws_zi_route <- tryCatch(posterior::as_draws_matrix(fit$draws("zi_route")),
                             error = function(e) NULL)
  draws_zi_obs <- tryCatch(posterior::as_draws_matrix(fit$draws("zi_obs")),
                           error = function(e) NULL)

  idx <- sample.int(nrow(draws_E), size = min(n_ppc_draws, nrow(draws_E)))
  y <- stan_data$count

  stat_fns <- list(
    mean = function(x) mean(x),
    sd = function(x) stats::sd(x),
    prop_zero = function(x) mean(x == 0),
    max = function(x) max(x)
  )

  obs_stats <- sapply(stat_fns, function(f) f(y))
  rep_stats <- matrix(NA_real_, nrow = length(idx), ncol = length(stat_fns),
                      dimnames = list(NULL, names(stat_fns)))

  for (j in seq_along(idx)) {
    d <- idx[j]
    mu <- exp(draws_E[d, ])
    phi <- as.numeric(draws_phi[d, 1])
    y_rep <- rnbinom(length(y), size = phi, mu = mu)

    if (!is.null(draws_zi_obs)) {
      structural_zero <- rbinom(length(y), size = 1, prob = draws_zi_obs[d, ]) == 1
      y_rep[structural_zero] <- 0L
    } else if (!is.null(draws_zi_route)) {
      zi_by_obs <- draws_zi_route[d, stan_data$route]
      structural_zero <- rbinom(length(y), size = 1, prob = zi_by_obs) == 1
      y_rep[structural_zero] <- 0L
    } else if (!is.null(draws_zi)) {
      zi <- as.numeric(draws_zi[d, 1])
      structural_zero <- rbinom(length(y), size = 1, prob = zi) == 1
      y_rep[structural_zero] <- 0L
    }

    rep_stats[j, ] <- sapply(stat_fns, function(f) f(y_rep))
  }

  bayes_p <- sapply(names(stat_fns), function(nm) {
    mean(rep_stats[, nm] >= obs_stats[nm])
  })

  ppc_long <- as.data.frame(rep_stats) %>%
    pivot_longer(everything(), names_to = "statistic", values_to = "value")
  obs_df <- tibble(statistic = names(obs_stats), value = as.numeric(obs_stats))

  p <- ggplot(ppc_long, aes(x = value)) +
    geom_histogram(bins = 40, fill = "#0072B2", alpha = 0.7) +
    geom_vline(data = obs_df, aes(xintercept = value),
               colour = "#D55E00", linewidth = 1) +
    facet_wrap(~ statistic, scales = "free") +
    labs(title = paste0("PPC: ", sp, " (", sp_code, ")"),
         subtitle = paste0("Model scenario: ", model_scenario),
         x = "Statistic value", y = "Replicate count") +
    theme_minimal()

  ggsave(file.path(diag_dir, paste0("ppc_", sp_f, "_", model_scenario, "_",
                                    land_cover, ".png")),
         p, width = 8, height = 6, dpi = 150)

  tibble(
    species = sp,
    species_code = sp_code,
    model_scenario = model_scenario,
    obs_mean = obs_stats["mean"],
    p_mean = bayes_p["mean"],
    obs_sd = obs_stats["sd"],
    p_sd = bayes_p["sd"],
    obs_prop_zero = obs_stats["prop_zero"],
    p_prop_zero = bayes_p["prop_zero"],
    obs_max = obs_stats["max"],
    p_max = bayes_p["max"]
  )
}

loo_summary <- function(fit, sp, sp_code, model_scenario) {
  if (!requireNamespace("loo", quietly = TRUE)) return(NULL)

  ll <- tryCatch(posterior::as_draws_matrix(fit$draws("log_lik")),
                 error = function(e) NULL)
  if (is.null(ll)) return(NULL)

  loo_fit <- loo::loo(ll)
  tibble(
    species = sp,
    species_code = sp_code,
    model_scenario = model_scenario,
    elpd_loo = loo_fit$estimates["elpd_loo", "Estimate"],
    p_loo = loo_fit$estimates["p_loo", "Estimate"],
    looic = loo_fit$estimates["looic", "Estimate"],
    pareto_k_bad = sum(loo_fit$diagnostics$pareto_k > 0.7, na.rm = TRUE)
  )
}

flag_ppc <- function(p) {
  ifelse(p < 0.05 | p > 0.95, "CHECK", "ok")
}

# Main ---------------------------------------------------------------------
manifest <- read.csv(manifest_file)

conv_rows <- list()
sampler_rows <- list()
ppc_rows <- list()
loo_rows <- list()
report_lines <- c(
  "##################################################",
  "# Scenario-aware Bayesian model diagnostics",
  paste("# Group:", land_cover, "| Years:", firstYear, "-", lastYear),
  paste("# Targets: R-hat <", rhat_thresh, "| ESS >", ess_thresh,
        "| divergences = 0 | E-BFMI >", ebfmi_thresh),
  "##################################################"
)

for (i in seq_len(nrow(manifest))) {
  sp <- manifest$Common.Name[i]
  sp_code <- manifest$Code[i]
  sp_f <- species_to_f(sp)
  model_scenario <- manifest$model_scenario[i]
  min_max_route_years <- manifest$min_max_route_years[i]
  base <- species_base(sp_f, model_scenario, min_max_route_years)

  fit_file <- file.path(fit_dir, paste0(base, "_stanfit.rds"))
  summary_file <- file.path(fit_dir, paste0(base, "_summ_fit.rds"))
  stan_data_file <- file.path(stan_data_dir, paste0(base, "_stan_data.RData"))

  report_lines <- c(
    report_lines,
    "",
    "==========================================================",
    paste0(" ", sp, " (", sp_code, ") | ", model_scenario),
    "=========================================================="
  )

  if (!file.exists(fit_file) && !file.exists(summary_file)) {
    report_lines <- c(report_lines, "  No saved fit found; skipping.")
    next
  }

  summ <- NULL
  if (file.exists(summary_file)) {
    summ <- readRDS(summary_file)
    cr <- convergence_summary(summ, sp, sp_code, model_scenario)
    if (!is.null(cr)) {
      conv_rows[[paste(sp_code, model_scenario, sep = "_")]] <- cr
      report_lines <- c(
        report_lines,
        "[1] Convergence",
        paste("    max R-hat        :", round(cr$max_rhat, 4),
              ifelse(cr$max_rhat <= rhat_thresh, "(ok)", "(CHECK)")),
        paste("    params R-hat >", rhat_thresh, ":", cr$n_rhat_gt_thresh),
        paste("    worst param      :", cr$worst_rhat_param),
        paste("    min bulk-ESS     :", round(cr$min_ess_bulk),
              ifelse(cr$min_ess_bulk >= ess_thresh, "(ok)", "(CHECK)")),
        paste("    min tail-ESS     :", round(cr$min_ess_tail),
              ifelse(cr$min_ess_tail >= ess_thresh, "(ok)", "(CHECK)"))
      )
    }
  }

  if (!file.exists(fit_file)) {
    report_lines <- c(report_lines, "[2/3/4] No full stanfit; sampler/PPC/LOO skipped.")
    next
  }

  fit <- tryCatch(readRDS(fit_file), error = function(e) NULL)
  if (is.null(fit)) {
    report_lines <- c(report_lines, "[2/3/4] stanfit could not be read; skipped.")
    next
  }
  if (!fit_output_files_exist(fit)) {
    report_lines <- c(report_lines, "[2/3/4] stanfit points to missing CmdStan CSV files; skipped.")
    next
  }

  sd_row <- sampler_diagnostics(fit, sp, sp_code, model_scenario)
  if (!is.null(sd_row)) {
    sampler_rows[[paste(sp_code, model_scenario, sep = "_")]] <- sd_row
    report_lines <- c(
      report_lines,
      "[2] Sampler",
      paste("    divergences      :", sd_row$num_divergent,
            ifelse(sd_row$num_divergent == 0, "(ok)", "(CHECK)")),
      paste("    max-treedepth    :", sd_row$num_max_treedepth),
      paste("    min E-BFMI       :", round(sd_row$min_ebfmi, 3),
            ifelse(sd_row$min_ebfmi >= ebfmi_thresh, "(ok)", "(CHECK)"))
    )
  }

  if (file.exists(stan_data_file)) {
    e <- new.env()
    load(stan_data_file, envir = e)
    ppc <- tryCatch(
      posterior_predictive_check(fit, e$stan_data, sp, sp_code, model_scenario, sp_f),
      error = function(err) {
        report_lines <<- c(report_lines, paste("[3] PPC failed:", conditionMessage(err)))
        NULL
      }
    )
    if (!is.null(ppc)) {
      ppc_rows[[paste(sp_code, model_scenario, sep = "_")]] <- ppc
      report_lines <- c(
        report_lines,
        "[3] PPC Bayesian p-values",
        paste("    mean      :", round(ppc$p_mean, 3), "(", flag_ppc(ppc$p_mean), ")"),
        paste("    sd        :", round(ppc$p_sd, 3), "(", flag_ppc(ppc$p_sd), ")"),
        paste("    prop zero :", round(ppc$p_prop_zero, 3), "(", flag_ppc(ppc$p_prop_zero), ")"),
        paste("    max       :", round(ppc$p_max, 3), "(", flag_ppc(ppc$p_max), ")")
      )
    }
  } else {
    report_lines <- c(report_lines, "[3] No stan_data file; PPC skipped.")
  }

  loo <- tryCatch(
    loo_summary(fit, sp, sp_code, model_scenario),
    error = function(err) {
      report_lines <<- c(report_lines, paste("[4] LOO failed:", conditionMessage(err)))
      NULL
    }
  )
  if (!is.null(loo)) {
    loo_rows[[paste(sp_code, model_scenario, sep = "_")]] <- loo
    report_lines <- c(
      report_lines,
      "[4] LOO",
      paste("    elpd_loo     :", round(loo$elpd_loo, 1)),
      paste("    pareto_k > .7:", loo$pareto_k_bad)
    )
  }

  rm(fit)
  gc()
}

writeLines(report_lines,
           file.path(diag_dir, paste0(land_cover, "_scenario_diagnostics_report.txt")))

conv_tbl <- bind_rows(conv_rows)
sampler_tbl <- bind_rows(sampler_rows)
ppc_tbl <- bind_rows(ppc_rows)
loo_tbl <- bind_rows(loo_rows)

if (nrow(conv_tbl) > 0) {
  write.csv(conv_tbl,
            file.path(diag_dir, paste0(land_cover, "_scenario_convergence_summary.csv")),
            row.names = FALSE)
}
if (nrow(sampler_tbl) > 0) {
  write.csv(sampler_tbl,
            file.path(diag_dir, paste0(land_cover, "_scenario_sampler_diagnostics.csv")),
            row.names = FALSE)
}
if (nrow(ppc_tbl) > 0) {
  write.csv(ppc_tbl,
            file.path(diag_dir, paste0(land_cover, "_scenario_ppc_summary.csv")),
            row.names = FALSE)
}
if (nrow(loo_tbl) > 0) {
  write.csv(loo_tbl,
            file.path(diag_dir, paste0(land_cover, "_scenario_loo_summary.csv")),
            row.names = FALSE)
}

combined <- conv_tbl %>%
  { if (ncol(.) == 0) tibble() else . }

if (nrow(combined) > 0 && nrow(sampler_tbl) > 0) {
  combined <- combined %>%
    left_join(sampler_tbl, by = c("species", "species_code", "model_scenario"))
}
if (nrow(combined) > 0 && nrow(ppc_tbl) > 0) {
  combined <- combined %>%
    left_join(ppc_tbl, by = c("species", "species_code", "model_scenario"))
}
if (nrow(combined) > 0 && nrow(loo_tbl) > 0) {
  combined <- combined %>%
    left_join(loo_tbl, by = c("species", "species_code", "model_scenario"))
}

if (nrow(combined) > 0) {
  write.csv(combined,
            file.path(diag_dir, paste0(land_cover, "_scenario_model_comparison_all.csv")),
            row.names = FALSE)
} else {
  message("No convergence rows were produced; combined comparison table skipped.")
}

cat("\n=== Scenario diagnostics complete ===\n")
cat("Outputs in:", diag_dir, "\n")
