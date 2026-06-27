# Test convergence / diagnostic improvements for the route-level habitat model.
#
# This script is intentionally separate from the production pipeline. It tests:
#   1. stricter route inclusion: min_max_route_years = 3
#   2. standardized habitat covariates
#   3. tighter scale priors in the lightweight NB test model
#   4. stable non-zero initial values
#   5. higher adapt_delta
#   6. a simple ZINB variant for species with excess-zero PPC failures
#
# Outputs are written under output/refit_tests/.

library(bbsBayes2)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(sf)
library(here)

here::i_am("code/test_convergence_refit_strategy.R")
source(here::here("functions", "neighbours_define_voronoi.R"))

# Settings -----------------------------------------------------------------
land_cover <- "grasslands"
legacy_land_cover <- "grassland"  # used only to read old diagnostics labels
firstYear <- 2010
lastYear <- 2024
year_range <- firstYear:lastYear
strat <- "bbs_usgs"
spatial_intercept <- TRUE

# Small, diagnostic-rich test panel:
#   BOSP = clean control, small and fast
#   CCLO = marginal convergence
#   UPSA = poor convergence, moderate size
#   BOBO = failed convergence, smaller than the very largest failed species
test_species_codes <- c("BOSP", "CCLO", "UPSA", "BOBO")

# For a very quick smoke test, set this to c("BOSP", "CCLO").
# For a stronger stress test, add "HOLA" or "EAKI", but those are slow.

min_max_route_years <- 3
fit_spatial_flag <- ifelse(spatial_intercept, 1L, 0L)

fit_models <- c("NB_test", "ZINB_test")

iter_warmup <- 1500
iter_sampling <- 1500
chains <- 4
parallel_chains <- 4
adapt_delta <- 0.99
max_treedepth <- 15
n_ppc_draws <- 250

# Set TRUE when you want to force fresh data prep. Fresh prep is needed to test
# min_max_route_years. If FALSE, the script will reuse a matching test cache.
rebuild_stan_data <- TRUE

out_dir <- here::here("output", "refit_tests")
data_cache_dir <- here::here("data", "stan_data", "refit_tests")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_cache_dir, recursive = TRUE, showWarnings = FALSE)

cov_csv <- here::here("data", paste0(land_cover, ".csv"))
mean_col <- paste0("mean_", land_cover)
slope_col <- paste0("slope_", land_cover)

nb_model_file <- here::here("models", "slope_habitat_route_NB_test.stan")
zinb_model_file <- here::here("models", "slope_habitat_route_ZINB_test.stan")

species_to_f <- function(sp) {
  gsub("'", "", gsub(" ", "_", sp, fixed = TRUE), fixed = TRUE)
}

standardize_safe <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) {
    return(as.numeric(x - mean(x, na.rm = TRUE)))
  }
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

sl <- function(y, x) {
  coefficients(stats::lm(y ~ x))[["x"]]
}

make_init <- function(stan_data, model_name) {
  force(stan_data)
  force(model_name)
  function() {
    init <- list(
      beta_raw = rnorm(stan_data$nroutes, 0, 0.02),
      BETA = 0,
      rho_beta_raw_hab = rnorm(stan_data$nroutes, 0, 0.02),
      rho_BETA_hab = 0,
      alpha_raw = rnorm(stan_data$nroutes, 0, 0.02),
      ALPHA = log(mean(stan_data$count + 0.5)),
      rho_alpha_raw_hab = rnorm(stan_data$nroutes, 0, 0.02),
      rho_ALPHA_hab = 0,
      eta = 0,
      obs_raw = rnorm(stan_data$nobservers, 0, 0.02),
      sdnoise = 0.5,
      sdobs = 0.1,
      sdbeta = 0.02,
      sdrho_beta_hab = 0.02,
      sdalpha = 0.5,
      sdrho_alpha_hab = 0.2
    )
    if (identical(model_name, "ZINB_test")) init$zi_logit <- -2
    init
  }
}

prepare_species_data <- function(species, species_bbs, species_f) {
  cache_file <- file.path(
    data_cache_dir,
    paste0(species_f, "_", land_cover, "_", firstYear, "_", lastYear,
           "_minmax", min_max_route_years, "_stan_data.RData")
  )

  if (!rebuild_stan_data && file.exists(cache_file)) {
    e <- new.env()
    load(cache_file, envir = e)
    return(as.list(e))
  }

  cat("\nPreparing data:", species, "\n")

  data_pkg <- bbsBayes2::stratify(by = strat, species = species_bbs,
                                  use_map = FALSE) %>%
    bbsBayes2::prepare_data(
      min_year = firstYear,
      max_year = lastYear,
      min_n_routes = 1,
      min_max_route_years = min_max_route_years
    )

  raw_data <- data_pkg[["raw_data"]] %>%
    filter(country_num == 840, state_num != 3)

  strata_map <- load_map(strat)

  route_map1 <- raw_data %>%
    select(route, strata_name, latitude, longitude) %>%
    distinct()

  route_map1 <- st_as_sf(route_map1, coords = c("longitude", "latitude"))
  st_crs(route_map1) <- 4326
  route_map1 <- st_transform(route_map1, crs = st_crs(strata_map))

  strata_map_buf <- strata_map %>%
    filter(strata_name %in% route_map1$strata_name) %>%
    summarise() %>%
    st_buffer(10000)

  realized_routes <- route_map1 %>%
    st_join(strata_map_buf, join = st_within, left = FALSE)

  new_data <- data.frame(
    strat_name = raw_data$strata_name,
    strat = raw_data$strata,
    route = raw_data$route,
    latitude = raw_data$latitude,
    longitude = raw_data$longitude,
    count = raw_data$count,
    year = raw_data$year_num,
    firstyr = raw_data$first_year,
    ObsN = raw_data$observer,
    r_year = raw_data$year
  ) %>%
    filter(route %in% realized_routes$route)

  cov_raw <- read.csv(cov_csv) %>%
    mutate(rt.uni = paste(StateNum, Route, sep = "-"))

  cov_full <- cov_raw %>%
    filter(!is.na(.data[[land_cover]])) %>%
    rename(route = rt.uni) %>%
    transmute(route, year, cov_value = .data[[land_cover]]) %>%
    filter(year %in% year_range)

  new_data <- new_data %>%
    inner_join(cov_full, by = c("route", "r_year" = "year"))

  if (nrow(new_data) == 0) {
    stop("No BBS count rows remain after route/covariate filtering for ", species)
  }

  strata_list <- data.frame(strata_name = unique(new_data$strat_name),
                            strat = unique(new_data$strat))
  realized_strata_map <- strata_map %>%
    filter(strata_name %in% strata_list$strata_name)

  new_data$routeF <- as.integer(factor(new_data$route))

  route_map <- unique(data.frame(
    route = new_data$route,
    routeF = new_data$routeF,
    strat = new_data$strat_name,
    latitude = new_data$latitude,
    longitude = new_data$longitude
  ))

  dups <- which(duplicated(route_map[, c("latitude", "longitude")]))
  while (length(dups) > 0) {
    route_map[dups, "latitude"] <- route_map[dups, "latitude"] + 0.01
    route_map[dups, "longitude"] <- route_map[dups, "longitude"] + 0.01
    dups <- which(duplicated(route_map[, c("latitude", "longitude")]))
  }

  route_map <- st_as_sf(route_map, coords = c("longitude", "latitude"))
  st_crs(route_map) <- 4326
  route_map <- st_transform(route_map, crs = st_crs(strata_map))

  car_stan_dat <- neighbours_define_voronoi(
    real_point_map = route_map,
    species = species,
    strat_indicator = "routeF",
    strata_map = realized_strata_map,
    concavity = 1
  )

  route_df <- route_map %>% sf::st_drop_geometry()

  slope_full <- cov_full %>%
    filter(route %in% route_df$route) %>%
    group_by(route) %>%
    summarise(
      !!mean_col := mean(cov_value, na.rm = TRUE),
      !!slope_col := sl(cov_value, year),
      .groups = "drop"
    ) %>%
    left_join(route_df, by = "route") %>%
    arrange(routeF) %>%
    filter(!is.na(routeF))

  stopifnot(nrow(slope_full) == max(new_data$routeF))
  stopifnot(identical(slope_full$routeF, seq_len(nrow(slope_full))))

  stan_data <- list(
    nroutes = max(new_data$routeF),
    ncounts = nrow(new_data),
    nyears = max(new_data$year),
    nobservers = max(as.integer(factor(new_data$ObsN))),
    count = as.integer(new_data$count),
    year = as.integer(new_data$year),
    route = as.integer(new_data$routeF),
    firstyr = as.integer(new_data$firstyr),
    observer = as.integer(factor(new_data$ObsN)),
    route_habitat = standardize_safe(slope_full[[mean_col]]),
    route_habitat_slope = standardize_safe(slope_full[[slope_col]]),
    fixedyear = as.integer(floor(stats::median(1:max(new_data$year)))),
    N_edges = car_stan_dat$N_edges,
    node1 = as.integer(car_stan_dat$node1),
    node2 = as.integer(car_stan_dat$node2),
    fit_spatial = fit_spatial_flag
  )

  if (car_stan_dat$N != stan_data$nroutes) {
    stop("Some routes are missing from adjacency matrix for ", species)
  }

  save(stan_data, new_data, route_map, realized_strata_map, car_stan_dat,
       slope_full, file = cache_file)

  list(
    stan_data = stan_data,
    new_data = new_data,
    route_map = route_map,
    realized_strata_map = realized_strata_map,
    car_stan_dat = car_stan_dat,
    slope_full = slope_full
  )
}

convergence_row <- function(fit, species, species_code, model_name) {
  s <- fit$summary() %>%
    filter(!is.na(rhat)) %>%
    filter(!str_detect(variable, "^(E|log_lik)\\["))

  tibble(
    species = species,
    species_code = species_code,
    model = model_name,
    n_parameters_checked = nrow(s),
    max_rhat = max(s$rhat, na.rm = TRUE),
    min_ess_bulk = min(s$ess_bulk, na.rm = TRUE),
    min_ess_tail = min(s$ess_tail, na.rm = TRUE),
    n_rhat_gt_1_01 = sum(s$rhat > 1.01, na.rm = TRUE),
    n_ess_lt_400 = sum(pmin(s$ess_bulk, s$ess_tail, na.rm = TRUE) < 400,
                       na.rm = TRUE),
    worst_rhat_param = s$variable[which.max(s$rhat)]
  )
}

sampler_row <- function(fit, species, species_code, model_name) {
  ds <- fit$diagnostic_summary(quiet = TRUE)
  tibble(
    species = species,
    species_code = species_code,
    model = model_name,
    num_divergent = sum(ds$num_divergent),
    num_max_treedepth = sum(ds$num_max_treedepth),
    min_ebfmi = min(ds$ebfmi, na.rm = TRUE)
  )
}

ppc_row <- function(fit, stan_data, species, species_code, model_name) {
  draws_E <- posterior::as_draws_matrix(fit$draws("E"))
  draws_phi <- posterior::as_draws_matrix(fit$draws("phi"))
  has_zi <- "zi" %in% variables(fit$draws())
  draws_zi <- if (has_zi) posterior::as_draws_matrix(fit$draws("zi")) else NULL

  idx <- sample.int(nrow(draws_E), min(n_ppc_draws, nrow(draws_E)))
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

    if (has_zi) {
      zi <- as.numeric(draws_zi[d, 1])
      structural_zero <- rbinom(length(y), size = 1, prob = zi) == 1
      y_rep[structural_zero] <- 0L
    }

    rep_stats[j, ] <- sapply(stat_fns, function(f) f(y_rep))
  }

  pvals <- sapply(names(stat_fns), function(nm) {
    mean(rep_stats[, nm] >= obs_stats[nm])
  })

  tibble(
    species = species,
    species_code = species_code,
    model = model_name,
    obs_mean = obs_stats["mean"],
    p_mean = pvals["mean"],
    obs_sd = obs_stats["sd"],
    p_sd = pvals["sd"],
    obs_prop_zero = obs_stats["prop_zero"],
    p_prop_zero = pvals["prop_zero"],
    obs_max = obs_stats["max"],
    p_max = pvals["max"]
  )
}

loo_row <- function(fit, species, species_code, model_name) {
  if (!requireNamespace("loo", quietly = TRUE)) return(NULL)
  ll <- posterior::as_draws_matrix(fit$draws("log_lik"))
  loo_fit <- loo::loo(ll)
  tibble(
    species = species,
    species_code = species_code,
    model = model_name,
    elpd_loo = loo_fit$estimates["elpd_loo", "Estimate"],
    p_loo = loo_fit$estimates["p_loo", "Estimate"],
    looic = loo_fit$estimates["looic", "Estimate"],
    pareto_k_bad = sum(loo_fit$diagnostics$pareto_k > 0.7, na.rm = TRUE)
  )
}

# Main ---------------------------------------------------------------------
spp_df <- read.csv(here::here("data", "spp_names_codes_group_aou.csv"))
target_spp <- spp_df %>%
  filter(Group == land_cover, Code %in% test_species_codes) %>%
  distinct(Common.Name, Code, .keep_all = TRUE) %>%
  mutate(test_order = match(Code, test_species_codes)) %>%
  arrange(test_order)

if (nrow(target_spp) == 0) {
  stop("No selected test species found for land_cover = ", land_cover)
}

cat("Test species:\n")
print(target_spp %>% select(Common.Name, Code, bbs_english))

mod_nb <- cmdstan_model(nb_model_file)
mod_zinb <- cmdstan_model(zinb_model_file)
models <- list(NB_test = mod_nb, ZINB_test = mod_zinb)

conv_rows <- list()
sampler_rows <- list()
ppc_rows <- list()
loo_rows <- list()

for (i in seq_len(nrow(target_spp))) {
  sp <- target_spp$Common.Name[i]
  sp_code <- target_spp$Code[i]
  sp_bbs <- target_spp$bbs_english[i]
  sp_f <- species_to_f(sp)

  prepared <- prepare_species_data(sp, sp_bbs, sp_f)
  stan_data <- prepared$stan_data

  cat("\n============================================================\n")
  cat(sp, "(", sp_code, "):", stan_data$nroutes, "routes,",
      stan_data$ncounts, "counts\n")
  cat("============================================================\n")

  for (model_name in fit_models) {
    cat("\nFitting", model_name, "for", sp, "\n")
    fit <- models[[model_name]]$sample(
      data = stan_data,
      seed = 1000 + i,
      chains = chains,
      parallel_chains = parallel_chains,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth,
      init = make_init(stan_data, model_name),
      refresh = 250,
      show_exceptions = TRUE,
      save_cmdstan_config = TRUE
    )

    base <- paste0(sp_f, "_", land_cover, "_", model_name,
                   "_minmax", min_max_route_years)
    saveRDS(fit, file.path(out_dir, paste0(base, "_stanfit.rds")))
    saveRDS(fit$summary(), file.path(out_dir, paste0(base, "_summary.rds")))

    conv_rows[[paste(sp_code, model_name, sep = "_")]] <-
      convergence_row(fit, sp, sp_code, model_name)
    sampler_rows[[paste(sp_code, model_name, sep = "_")]] <-
      sampler_row(fit, sp, sp_code, model_name)
    ppc_rows[[paste(sp_code, model_name, sep = "_")]] <-
      ppc_row(fit, stan_data, sp, sp_code, model_name)
    loo_rows[[paste(sp_code, model_name, sep = "_")]] <-
      loo_row(fit, sp, sp_code, model_name)

    rm(fit)
    gc()
  }

  conv_tbl <- bind_rows(conv_rows)
  sampler_tbl <- bind_rows(sampler_rows)
  ppc_tbl <- bind_rows(ppc_rows)
  loo_tbl <- bind_rows(loo_rows)

  write.csv(conv_tbl, file.path(out_dir, "test_convergence_summary.csv"),
            row.names = FALSE)
  write.csv(sampler_tbl, file.path(out_dir, "test_sampler_diagnostics.csv"),
            row.names = FALSE)
  write.csv(ppc_tbl, file.path(out_dir, "test_ppc_summary.csv"),
            row.names = FALSE)
  if (nrow(loo_tbl) > 0) {
    write.csv(loo_tbl, file.path(out_dir, "test_loo_summary.csv"),
              row.names = FALSE)
  }
}

comparison <- bind_rows(conv_rows) %>%
  left_join(bind_rows(sampler_rows),
            by = c("species", "species_code", "model")) %>%
  left_join(bind_rows(ppc_rows),
            by = c("species", "species_code", "model"))

loo_tbl <- bind_rows(loo_rows)
if (nrow(loo_tbl) > 0) {
  comparison <- comparison %>%
    left_join(loo_tbl, by = c("species", "species_code", "model"))
}

write.csv(comparison, file.path(out_dir, "test_model_comparison_all.csv"),
          row.names = FALSE)

cat("\n=== Refit strategy test complete ===\n")
cat("Outputs:", out_dir, "\n")
print(comparison %>%
        select(species_code, model, max_rhat, min_ess_bulk, num_divergent,
               min_ebfmi, p_prop_zero, p_sd, p_max, any_of("elpd_loo")))
