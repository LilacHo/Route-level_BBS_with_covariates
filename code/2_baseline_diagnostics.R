# 2_baseline_diagnostics.R
#
# Reliability / fit diagnostics for the per-species Bayesian hierarchical
# route-level slope models fit in 1_baseline_NB_fit.R. These baseline
# diagnostics feed the model-scenario assignment in 3_assign_model_scenarios.R.
#
# These are NOT "validation" against held-out truth; they are the standard
# checks used to decide whether the posterior can be trusted and whether the
# negative-binomial model adequately reproduces the observed counts. The
# modelling framework follows Smith et al. 2024 (https://ace-eco.org/vol19/iss2/art23/).
#
# Three families of checks are produced for every species in the target group:
#
#   1. MCMC convergence        : R-hat, bulk-ESS, tail-ESS (posterior pkg).
#                                R-hat < 1.01 and ESS > 400 are the targets
#                                (Vehtari et al. 2021).
#   2. HMC/NUTS sampler health : divergent transitions (want 0), max-treedepth
#                                saturation, and E-BFMI (want > 0.2).
#   3. Posterior predictive    : simulate replicate counts from the fitted
#      checks (PPC)              negative binomial and compare summary stats
#                                (mean, sd, proportion of zeros, max) to the
#                                observed counts via Bayesian p-values.
#
# Inputs  (per species, written by 1_baseline_NB_fit.R):
#   output/<species>_<land_cover>_<firstYear>_<lastYear>_stanfit.rds  (CmdStanMCMC)
#   output/<species>_<land_cover>_<firstYear>_<lastYear>_summ_fit.rds (draws summary)
#   data/stan_data/<species>_<land_cover>_<firstYear>_<lastYear>_stan_data.RData
#
# Outputs:
#   output/model_diagnostics/<land_cover>_convergence_summary.csv
#   output/model_diagnostics/<land_cover>_sampler_diagnostics.csv
#   output/model_diagnostics/<land_cover>_posterior_predictive_checks.csv
#   output/model_diagnostics/<land_cover>_diagnostics_report.txt
#   output/model_diagnostics/ppc_<species>_<land_cover>.png

library(here)
library(tidyverse)
library(posterior)
library(cmdstanr)

here::i_am("code/2_baseline_diagnostics.R")

# Settings -----------------------------------------------------------------
land_cover <- "grasslands"   # must match the group fit in script 1
firstYear  <- 2010
lastYear   <- 2024

# Diagnostic thresholds (Vehtari et al. 2021 / Stan recommendations)
rhat_thresh   <- 1.01
ess_thresh    <- 400
ebfmi_thresh  <- 0.2
n_ppc_draws   <- 500   # posterior draws subsampled for the predictive check

fit_dir   <- here::here("output")
diag_dir  <- here::here("output", "model_diagnostics")
stan_data_dir <- here::here("data", "stan_data")
if (!dir.exists(diag_dir)) dir.create(diag_dir, recursive = TRUE)

# Species list for the target group ----------------------------------------
spp_df <- read.csv(here::here("data", "spp_names_codes_group_aou.csv"))

species_to_f <- function(sp) {
  gsub("'", "", gsub(" ", "_", sp, fixed = TRUE), fixed = TRUE)
}

target_spp <- spp_df %>%
  filter(Group == land_cover) %>%
  distinct(Common.Name, Code, .keep_all = TRUE)

if (nrow(target_spp) == 0) {
  stop("No species found for land_cover = '", land_cover, "' in spp_names_codes_group_aou.csv")
}


# ==========================================================================
# Functions ####
# ==========================================================================

# Build the per-species output basename used by script 1.
species_base <- function(sp_f) {
  paste0(sp_f, "_", land_cover, "_", firstYear, "_", lastYear)
}

# ---- 1. Convergence summary from the saved draws summary ------------------
# The *_summ_fit.rds object is the output of CmdStanMCMC$summary(), which is a
# posterior::summarise_draws() data frame and already carries rhat/ess columns.
convergence_summary <- function(summ, sp, sp_code) {
  needed <- c("rhat", "ess_bulk", "ess_tail")
  if (!all(needed %in% names(summ))) {
    # Older summaries may lack these; nothing to report.
    return(NULL)
  }

  s <- summ %>% filter(!is.na(rhat))
  n_par <- nrow(s)

  tibble(
    species        = sp,
    species_code   = sp_code,
    n_parameters   = n_par,
    max_rhat       = max(s$rhat, na.rm = TRUE),
    n_rhat_gt_thresh = sum(s$rhat > rhat_thresh, na.rm = TRUE),
    prop_rhat_ok   = mean(s$rhat <= rhat_thresh, na.rm = TRUE),
    min_ess_bulk   = min(s$ess_bulk, na.rm = TRUE),
    min_ess_tail   = min(s$ess_tail, na.rm = TRUE),
    n_ess_lt_thresh = sum(pmin(s$ess_bulk, s$ess_tail, na.rm = TRUE) < ess_thresh,
                          na.rm = TRUE),
    # name of the worst-mixing parameter (largest rhat) for quick triage
    worst_rhat_param = s$variable[which.max(s$rhat)]
  )
}

# ---- 2. Sampler (HMC/NUTS) diagnostics from the CmdStanMCMC object --------
sampler_diagnostics <- function(stanfit, sp, sp_code) {
  # diagnostic_summary() returns counts of divergences / treedepth hits and
  # the per-chain E-BFMI; wrapped in tryCatch because it is only available for
  # MCMC fits.
  ds <- tryCatch(stanfit$diagnostic_summary(quiet = TRUE),
                 error = function(e) NULL)
  if (is.null(ds)) return(NULL)

  n_div   <- sum(ds$num_divergent)
  n_tree  <- sum(ds$num_max_treedepth)
  ebfmi   <- ds$ebfmi

  # Total post-warmup draws across chains (for context on divergence rate)
  n_draws <- tryCatch(posterior::ndraws(stanfit$draws()), error = function(e) NA_integer_)

  tibble(
    species          = sp,
    species_code     = sp_code,
    num_divergent    = n_div,
    divergence_rate  = if (!is.na(n_draws) && n_draws > 0) n_div / n_draws else NA_real_,
    num_max_treedepth = n_tree,
    min_ebfmi        = min(ebfmi, na.rm = TRUE),
    n_chains_low_ebfmi = sum(ebfmi < ebfmi_thresh, na.rm = TRUE)
  )
}

# ---- 3. Posterior predictive checks for the negative-binomial counts ------
# The model's likelihood is count ~ neg_binomial_2_log(E, phi). We draw a
# subset of posterior samples of (E, phi), simulate replicate count vectors,
# and compare four discrepancy statistics against the observed counts. The
# Bayesian p-value is the proportion of replicate stats exceeding the observed
# one; values near 0 or 1 indicate the model fails to reproduce that feature.
posterior_predictive_check <- function(stanfit, observed_counts, sp, sp_code, sp_f) {
  draws_E <- tryCatch(
    posterior::as_draws_matrix(stanfit$draws("E")),
    error = function(e) NULL
  )
  draws_phi <- tryCatch(
    posterior::as_draws_matrix(stanfit$draws("phi")),
    error = function(e) NULL
  )
  if (is.null(draws_E) || is.null(draws_phi)) return(NULL)

  # Zero-inflation probability (theta) is only present in the ZINB model. When
  # it exists we generate structural zeros before the NB draw so the predictive
  # check matches the fitted likelihood; when absent theta = 0 (plain NB).
  draws_theta <- tryCatch(
    posterior::as_draws_matrix(stanfit$draws("theta")),
    error = function(e) NULL
  )

  total_draws <- nrow(draws_E)
  idx <- sample.int(total_draws, size = min(n_ppc_draws, total_draws))
  draws_E   <- draws_E[idx, , drop = FALSE]
  draws_phi <- as.numeric(draws_phi[idx, 1])
  draws_theta <- if (is.null(draws_theta)) rep(0, length(idx)) else as.numeric(draws_theta[idx, 1])

  n_sim <- length(idx)
  ncounts <- length(observed_counts)

  # Discrepancy statistics computed on each replicate dataset.
  stat_fns <- list(
    mean      = function(y) mean(y),
    sd        = function(y) sd(y),
    prop_zero = function(y) mean(y == 0),
    max       = function(y) max(y)
  )
  obs_stats <- sapply(stat_fns, function(f) f(observed_counts))
  rep_stats <- matrix(NA_real_, nrow = n_sim, ncol = length(stat_fns),
                      dimnames = list(NULL, names(stat_fns)))

  for (i in seq_len(n_sim)) {
    mu  <- exp(draws_E[i, ])               # expected count on natural scale
    # rnbinom 'size' = phi (dispersion), 'mu' parameterisation matches neg_binomial_2
    y_rep <- rnbinom(ncounts, size = draws_phi[i], mu = mu)
    # Apply zero-inflation: each observation is a structural zero with prob theta
    if (draws_theta[i] > 0) {
      structural_zero <- rbinom(ncounts, size = 1, prob = draws_theta[i])
      y_rep[structural_zero == 1] <- 0
    }
    rep_stats[i, ] <- sapply(stat_fns, function(f) f(y_rep))
  }

  # Bayesian p-value: P(T(y_rep) >= T(y_obs))
  bayes_p <- sapply(names(stat_fns), function(nm) {
    mean(rep_stats[, nm] >= obs_stats[nm])
  })

  # Density-overlay style plot of one statistic set: observed vs replicate
  ppc_long <- as.data.frame(rep_stats) %>%
    pivot_longer(everything(), names_to = "statistic", values_to = "value")
  obs_df <- tibble(statistic = names(obs_stats), value = as.numeric(obs_stats))

  p <- ggplot(ppc_long, aes(x = value)) +
    geom_histogram(bins = 40, fill = "#0072B2", alpha = 0.7) +
    geom_vline(data = obs_df, aes(xintercept = value),
               colour = "#D55E00", linewidth = 1) +
    facet_wrap(~ statistic, scales = "free") +
    labs(title = paste0("Posterior predictive check: ", sp, " (", sp_code, ")"),
         subtitle = "Histogram = replicate datasets; red line = observed",
         x = "Statistic value", y = "Count of replicate datasets") +
    theme_minimal()
  ggsave(file.path(diag_dir, paste0("ppc_", sp_f, "_", land_cover, ".png")),
         p, width = 8, height = 6, dpi = 150)

  tibble(
    species      = sp,
    species_code = sp_code,
    obs_mean      = obs_stats["mean"],      p_mean      = bayes_p["mean"],
    obs_sd        = obs_stats["sd"],        p_sd        = bayes_p["sd"],
    obs_prop_zero = obs_stats["prop_zero"], p_prop_zero = bayes_p["prop_zero"],
    obs_max       = obs_stats["max"],       p_max       = bayes_p["max"]
  )
}

# Flag p-values that indicate poor fit (extreme in either tail).
flag_ppc <- function(p) ifelse(p < 0.05 | p > 0.95, "CHECK", "ok")


# ==========================================================================
# Main loop ####
# ==========================================================================

conv_rows   <- list()
samp_rows   <- list()
ppc_rows    <- list()

report_file <- file.path(diag_dir, paste0(land_cover, "_diagnostics_report.txt"))
sink(report_file)

cat("##################################################\n")
cat("# Bayesian hierarchical model diagnostics\n")
cat("# Group:", land_cover, "| Years:", firstYear, "-", lastYear, "\n")
cat("# Targets: R-hat <", rhat_thresh, "| ESS >", ess_thresh,
    "| divergences = 0 | E-BFMI >", ebfmi_thresh, "\n")
cat("##################################################\n")

for (i in seq_len(nrow(target_spp))) {
  sp      <- target_spp$Common.Name[i]
  sp_f    <- species_to_f(sp)
  sp_code <- target_spp$Code[i]
  base    <- species_base(sp_f)

  summ_file    <- file.path(fit_dir, paste0(base, "_summ_fit.rds"))
  stanfit_file <- file.path(fit_dir, paste0(base, "_stanfit.rds"))

  cat("\n==========================================================\n")
  cat(" ", sp, " (", sp_code, ")\n")
  cat("==========================================================\n")

  if (!file.exists(summ_file) && !file.exists(stanfit_file)) {
    cat("  No saved fit found — skipping.\n")
    next
  }

  # ---- Convergence (from summary) ----
  if (file.exists(summ_file)) {
    summ <- readRDS(summ_file)
    cr <- convergence_summary(summ, sp, sp_code)
    if (!is.null(cr)) {
      conv_rows[[sp]] <- cr
      cat("\n[1] Convergence\n")
      cat("    max R-hat        :", round(cr$max_rhat, 4),
          ifelse(cr$max_rhat <= rhat_thresh, " (ok)", " (CHECK)"), "\n")
      cat("    params R-hat >",  rhat_thresh, ":", cr$n_rhat_gt_thresh, "\n")
      cat("    worst param      :", cr$worst_rhat_param, "\n")
      cat("    min bulk-ESS     :", round(cr$min_ess_bulk),
          ifelse(cr$min_ess_bulk >= ess_thresh, " (ok)", " (CHECK)"), "\n")
      cat("    min tail-ESS     :", round(cr$min_ess_tail),
          ifelse(cr$min_ess_tail >= ess_thresh, " (ok)", " (CHECK)"), "\n")
    } else {
      cat("\n[1] Convergence: summary lacks rhat/ess columns — skipped.\n")
    }
  } else {
    cat("\n[1] Convergence: no summary file — skipped.\n")
  }

  # ---- Sampler diagnostics + PPC (need the full stanfit) ----
  if (file.exists(stanfit_file)) {
    stanfit <- tryCatch(readRDS(stanfit_file), error = function(e) NULL)

    if (is.null(stanfit)) {
      cat("\n[2] Sampler diagnostics: stanfit could not be read — skipped.\n")
    } else {
      sd_row <- sampler_diagnostics(stanfit, sp, sp_code)
      if (!is.null(sd_row)) {
        samp_rows[[sp]] <- sd_row
        cat("\n[2] Sampler (HMC/NUTS)\n")
        cat("    divergences      :", sd_row$num_divergent,
            ifelse(sd_row$num_divergent == 0, " (ok)", " (CHECK)"), "\n")
        cat("    max-treedepth    :", sd_row$num_max_treedepth, "\n")
        cat("    min E-BFMI       :", round(sd_row$min_ebfmi, 3),
            ifelse(sd_row$min_ebfmi >= ebfmi_thresh, " (ok)", " (CHECK)"), "\n")
      }

      # PPC needs observed counts from the saved stan_data
      stan_data_file <- file.path(stan_data_dir, paste0(base, "_stan_data.RData"))
      if (file.exists(stan_data_file)) {
        e <- new.env()
        load(stan_data_file, envir = e)
        observed_counts <- e$stan_data$count
        pr <- tryCatch(
          posterior_predictive_check(stanfit, observed_counts, sp, sp_code, sp_f),
          error = function(err) {
            cat("\n[3] PPC failed:", conditionMessage(err), "\n"); NULL
          }
        )
        if (!is.null(pr)) {
          ppc_rows[[sp]] <- pr
          cat("\n[3] Posterior predictive checks (Bayesian p-values)\n")
          cat("    mean      : p =", round(pr$p_mean, 3),      "(", flag_ppc(pr$p_mean), ")\n")
          cat("    sd        : p =", round(pr$p_sd, 3),        "(", flag_ppc(pr$p_sd), ")\n")
          cat("    prop zero : p =", round(pr$p_prop_zero, 3), "(", flag_ppc(pr$p_prop_zero), ")\n")
          cat("    max       : p =", round(pr$p_max, 3),       "(", flag_ppc(pr$p_max), ")\n")
        }
      } else {
        cat("\n[3] PPC: no stan_data file — skipped.\n")
      }
    }
  } else {
    cat("\n[2/3] No stanfit file — sampler diagnostics and PPC skipped.\n")
  }
}

sink()
cat("Saved per-species diagnostics report to:", report_file, "\n")

# ---- Write combined CSV tables -------------------------------------------
if (length(conv_rows) > 0) {
  conv_tbl <- bind_rows(conv_rows)
  write.csv(conv_tbl,
            file.path(diag_dir, paste0(land_cover, "_convergence_summary.csv")),
            row.names = FALSE)
  cat("Saved convergence summary (", nrow(conv_tbl), "species )\n")
}

if (length(samp_rows) > 0) {
  samp_tbl <- bind_rows(samp_rows)
  write.csv(samp_tbl,
            file.path(diag_dir, paste0(land_cover, "_sampler_diagnostics.csv")),
            row.names = FALSE)
  cat("Saved sampler diagnostics (", nrow(samp_tbl), "species )\n")
}

if (length(ppc_rows) > 0) {
  ppc_tbl <- bind_rows(ppc_rows)
  write.csv(ppc_tbl,
            file.path(diag_dir, paste0(land_cover, "_posterior_predictive_checks.csv")),
            row.names = FALSE)
  cat("Saved posterior predictive checks (", nrow(ppc_tbl), "species )\n")
}

cat("\n=== Model diagnostics complete ===\n")
cat("Outputs in:", diag_dir, "\n")
