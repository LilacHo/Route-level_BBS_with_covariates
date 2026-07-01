# Fit route-level BBS habitat models and calculate route-level CH summaries.
#
# This scenario-aware version can assign different species to different Stan
# models. It keeps the same core output columns used downstream:
#   species, species_code, land_cover, route, routeF, latitude, longitude,
#   model_scenario, CH, CH_no_habitat, CH_dif
#
# The default scenario is the improved negative-binomial test model. Species
# with known excess-zero problems can be assigned to route-level or
# route-plus-abundance zero-inflated models in the "Species model scenarios"
# block below.

library(bbsBayes2)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(sf)
library(here)

here::i_am("code/1_CH_no_habitat_routes.R")
source(here::here("functions", "neighbours_define_voronoi.R"))

# Settings ----------------------------------------------------------------
land_cover <- "grasslands"
firstYear <- 2010
lastYear <- 2024
year_range <- firstYear:lastYear
dt <- lastYear - firstYear

strat <- "bbs_usgs"
spatial_intercept <- TRUE
fit_spatial_flag <- ifelse(spatial_intercept, 1L, 0L)

# If FALSE, existing scenario fit files are reused when possible.
refit_existing_models <- FALSE
rebuild_stan_data <- FALSE

# If TRUE, writes the selected scenario route summaries to the legacy
# output/species_routes folder used by steps 2 and 3.
write_primary_route_csv <- TRUE

# Optional quick-run subset. Wrapper scripts can define only_species_codes
# before sourcing this file; otherwise all species in land_cover are used.
if (!exists("only_species_codes")) {
  only_species_codes <- character(0)
}

# Scenario assignments produced by code/1a_prepare_model_scenarios.R. If this
# file exists, it overrides the fallback vectors below.
scenario_assignment_file <- here::here(
  "output", "model_scenario_assignments",
  paste0(land_cover, "_model_scenario_assignments.csv")
)

cov_csv <- here::here("data", paste0(land_cover, ".csv"))
mean_col <- paste0("mean_", land_cover)
slope_col <- paste0("slope_", land_cover)

fit_dir <- here::here("output", "model_fits_by_scenario")
route_out_dir <- here::here("output", "species_routes_by_scenario")
legacy_route_out_dir <- here::here("output", "species_routes")
stan_data_dir <- here::here("data", "stan_data_by_scenario")
cmdstan_csv_dir <- here::here("output", "cmdstan_csv_by_scenario")
dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(route_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(legacy_route_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stan_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cmdstan_csv_dir, recursive = TRUE, showWarnings = FALSE)

# Model catalog ------------------------------------------------------------
model_catalog <- tibble::tribble(
  ~model_scenario,       ~model_file,                                                   ~min_max_route_years, ~iter_warmup, ~iter_sampling, ~chains, ~adapt_delta, ~max_treedepth,
  "NB_test",             here::here("models", "slope_habitat_route_NB_test.stan"),       3L,                  2000L,        2000L,          4L,      0.99,         15L,
  "ZINB_route_test",     here::here("models", "slope_habitat_route_ZINB_route_test.stan"), 3L,                2000L,        2000L,          4L,      0.99,         15L,
  "ZINB_route_mu_test",  here::here("models", "slope_habitat_route_ZINB_route_mu_test.stan"), 3L,             3000L,        3000L,          4L,      0.99,         15L
)

# Fallback species model scenarios -----------------------------------------
# These are used only if scenario_assignment_file does not exist.
route_zinb_species <- c(
  "CCLO"
)

route_mu_zinb_species <- c(
  "BOBO", "UPSA",
  "CASP", "CCSP", "DICK", "EAKI", "EAME", "GRSP", "HOLA", "LARB",
  "LOSH", "NOBO", "RNEP", "SAVS", "SEWR", "STFL", "VESP", "WEKI", "WEME"
)

scenario_for_species <- function(species_code) {
  if (species_code %in% route_mu_zinb_species) return("ZINB_route_mu_test")
  if (species_code %in% route_zinb_species) return("ZINB_route_test")
  "NB_test"
}

load_scenario_assignments <- function() {
  if (!file.exists(scenario_assignment_file)) return(NULL)

  assignments <- read.csv(scenario_assignment_file)
  required <- c("Code", "model_scenario")
  if (!all(required %in% names(assignments))) {
    stop("Scenario assignment file must contain columns: ",
         paste(required, collapse = ", "))
  }
  assignments %>%
    select(Code, model_scenario,
           any_of(c("assignment_reason", "diagnostic_flags"))) %>%
    distinct(Code, .keep_all = TRUE)
}

# Helpers ------------------------------------------------------------------
species_to_f <- function(sp) {
  gsub("'", "", gsub(" ", "_", sp, fixed = TRUE), fixed = TRUE)
}

standardize_safe <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(as.numeric(x - mean(x, na.rm = TRUE)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

sl <- function(y, x) {
  coefficients(stats::lm(y ~ x))[["x"]]
}

pct_cum <- function(b, dt) {
  100 * (exp(b * dt) - 1)
}

make_init <- function(stan_data, model_scenario) {
  force(stan_data)
  force(model_scenario)
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

    if (identical(model_scenario, "ZINB_route_test")) {
      init$zi_intercept <- -2
      init$zi_route_raw <- rnorm(stan_data$nroutes, 0, 0.02)
      init$sd_zi_route <- 0.25
    }
    if (identical(model_scenario, "ZINB_route_mu_test")) {
      init$zi_intercept <- -2
      init$zi_log_mu <- -0.5
      init$zi_route_raw <- rnorm(stan_data$nroutes, 0, 0.02)
      init$sd_zi_route <- 0.25
    }
    init
  }
}

fit_output_files_exist <- function(fit) {
  files <- tryCatch(
    fit$output_files(include_failed = FALSE),
    error = function(e) character()
  )
  length(files) > 0 && all(file.exists(files))
}

save_fit_with_csv <- function(fit, fit_file, summary_file, csv_dir) {
  dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
  fit$save_output_files(dir = csv_dir)
  saveRDS(fit, fit_file)
  saveRDS(fit$summary(), summary_file)
}

prepare_species_data <- function(species, species_bbs, species_f, model_scenario,
                                 min_max_route_years) {
  cache_file <- file.path(
    stan_data_dir,
    paste0(species_f, "_", land_cover, "_", model_scenario, "_",
           firstYear, "_", lastYear, "_minmax", min_max_route_years,
           "_stan_data.RData")
  )

  if (!rebuild_stan_data && file.exists(cache_file)) {
    e <- new.env()
    load(cache_file, envir = e)
    return(as.list(e))
  }

  cat("  Preparing data with min_max_route_years =", min_max_route_years, "\n")

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
    stop("No BBS count rows remain after route/covariate filtering")
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

  observer <- as.integer(factor(new_data$ObsN))
  stan_data <- list(
    nroutes = max(new_data$routeF),
    ncounts = nrow(new_data),
    nyears = max(new_data$year),
    nobservers = max(observer),
    count = as.integer(new_data$count),
    year = as.integer(new_data$year),
    route = as.integer(new_data$routeF),
    firstyr = as.integer(new_data$firstyr),
    observer = observer,
    route_habitat = standardize_safe(slope_full[[mean_col]]),
    route_habitat_slope = standardize_safe(slope_full[[slope_col]]),
    fixedyear = as.integer(floor(stats::median(1:max(new_data$year)))),
    N_edges = car_stan_dat$N_edges,
    node1 = as.integer(car_stan_dat$node1),
    node2 = as.integer(car_stan_dat$node2),
    fit_spatial = fit_spatial_flag
  )

  if (car_stan_dat$N != stan_data$nroutes) {
    stop("Some routes are missing from adjacency matrix")
  }

  dist_matrix_km <- dist_matrix(route_map, strat_indicator = "routeF")

  save(stan_data, new_data, route_map, realized_strata_map, car_stan_dat,
       dist_matrix_km, cov_full, slope_full, land_cover, model_scenario,
       file = cache_file)

  list(
    stan_data = stan_data,
    new_data = new_data,
    route_map = route_map,
    realized_strata_map = realized_strata_map,
    car_stan_dat = car_stan_dat,
    dist_matrix_km = dist_matrix_km,
    cov_full = cov_full,
    slope_full = slope_full,
    land_cover = land_cover,
    model_scenario = model_scenario
  )
}

route_summary_from_fit <- function(summ, route_map, species, species_code,
                                   model_scenario) {
  beta_df <- summ %>%
    filter(str_detect(variable, "^beta\\[")) %>%
    transmute(routeF = as.integer(str_extract(variable, "\\d+")),
              beta = mean)

  beta_resid_df <- summ %>%
    filter(str_detect(variable, "^beta_resid\\[")) %>%
    transmute(routeF = as.integer(str_extract(variable, "\\d+")),
              beta_resid = mean)

  route_info <- route_map %>%
    st_transform(4326) %>%
    mutate(longitude = st_coordinates(.)[, 1],
           latitude = st_coordinates(.)[, 2]) %>%
    st_drop_geometry() %>%
    select(route, routeF, latitude, longitude)

  beta_df %>%
    left_join(beta_resid_df, by = "routeF") %>%
    mutate(
      CH = pct_cum(beta, dt),
      CH_no_habitat = pct_cum(beta_resid, dt),
      CH_dif = CH - CH_no_habitat
    ) %>%
    left_join(route_info, by = "routeF") %>%
    mutate(
      species = species,
      species_code = species_code,
      land_cover = land_cover,
      model_scenario = model_scenario
    ) %>%
    select(species, species_code, land_cover, model_scenario, route, routeF,
           latitude, longitude, CH, CH_no_habitat, CH_dif)
}

# Main ---------------------------------------------------------------------
spp_df <- read.csv(here::here("data", "spp_names_codes_group_aou.csv"))
scenario_assignments <- load_scenario_assignments()

target_spp <- spp_df %>%
  filter(Group == land_cover) %>%
  distinct(Common.Name, Code, .keep_all = TRUE)

if (!is.null(scenario_assignments)) {
  cat("Using model scenarios from:", scenario_assignment_file, "\n")
  target_spp <- target_spp %>%
    left_join(scenario_assignments, by = "Code") %>%
    mutate(model_scenario = ifelse(is.na(model_scenario),
                                   vapply(Code, scenario_for_species, character(1)),
                                   model_scenario))
} else {
  cat("No scenario assignment file found; using fallback scenario rules.\n")
  target_spp <- target_spp %>%
    mutate(model_scenario = vapply(Code, scenario_for_species, character(1)))
}

if (length(only_species_codes) > 0) {
  target_spp <- target_spp %>% filter(Code %in% only_species_codes)
}

target_spp <- target_spp %>%
  left_join(model_catalog, by = "model_scenario")

if (any(is.na(target_spp$model_file))) {
  stop("At least one species was assigned to an unknown model_scenario")
}

cat(land_cover, "species to fit (n =", nrow(target_spp), "):\n")
print(target_spp %>% select(Common.Name, Code, bbs_english, model_scenario))

compiled_models <- list()
results_list <- list()

for (i in seq_len(nrow(target_spp))) {
  sp <- target_spp$Common.Name[i]
  sp_f <- species_to_f(sp)
  sp_code <- target_spp$Code[i]
  sp_bbs <- target_spp$bbs_english[i]
  model_scenario <- target_spp$model_scenario[i]
  model_file <- target_spp$model_file[i]
  min_max_route_years <- target_spp$min_max_route_years[i]

  cat("\n[", i, "/", nrow(target_spp), "] ", sp, " (", sp_code, ") -> ",
      model_scenario, "\n", sep = "")

  base <- paste0(sp_f, "_", land_cover, "_", model_scenario, "_",
                 firstYear, "_", lastYear, "_minmax", min_max_route_years)
  fit_file <- file.path(fit_dir, paste0(base, "_stanfit.rds"))
  summary_file <- file.path(fit_dir, paste0(base, "_summ_fit.rds"))
  csv_dir <- file.path(cmdstan_csv_dir, base)
  route_csv <- file.path(route_out_dir,
                         paste0(land_cover, "_", model_scenario, "_", sp_f,
                                "_route_CH.csv"))
  legacy_route_csv <- file.path(legacy_route_out_dir,
                                paste0(land_cover, "_", sp_f, "_route_CH.csv"))

  prepared <- prepare_species_data(
    species = sp,
    species_bbs = sp_bbs,
    species_f = sp_f,
    model_scenario = model_scenario,
    min_max_route_years = min_max_route_years
  )
  stan_data <- prepared$stan_data
  route_map <- prepared$route_map

  if (file.exists(fit_file) && !refit_existing_models) {
    cat("  Loading existing scenario fit\n")
    stanfit <- readRDS(fit_file)
    if (!fit_output_files_exist(stanfit)) {
      stop("Existing fit points to missing CmdStan CSV files. Refit by setting ",
           "refit_existing_models <- TRUE or deleting: ", fit_file)
    }
    if (file.exists(summary_file)) {
      summ <- readRDS(summary_file)
    } else {
      summ <- stanfit$summary()
      saveRDS(summ, summary_file)
    }
  } else {
    if (is.null(compiled_models[[model_scenario]])) {
      compiled_models[[model_scenario]] <- cmdstan_model(model_file)
    }

    cat("  Fitting model\n")
    stanfit <- compiled_models[[model_scenario]]$sample(
      data = stan_data,
      seed = 1000 + i,
      chains = target_spp$chains[i],
      parallel_chains = target_spp$chains[i],
      iter_warmup = target_spp$iter_warmup[i],
      iter_sampling = target_spp$iter_sampling[i],
      adapt_delta = target_spp$adapt_delta[i],
      max_treedepth = target_spp$max_treedepth[i],
      init = make_init(stan_data, model_scenario),
      refresh = 400,
      show_exceptions = TRUE,
      save_cmdstan_config = TRUE
    )

    save_fit_with_csv(stanfit, fit_file, summary_file, csv_dir)
    summ <- readRDS(summary_file)
  }

  route_ch <- route_summary_from_fit(
    summ = summ,
    route_map = route_map,
    species = sp,
    species_code = sp_code,
    model_scenario = model_scenario
  )

  write.csv(route_ch, route_csv, row.names = FALSE)
  if (write_primary_route_csv) {
    write.csv(route_ch, legacy_route_csv, row.names = FALSE)
  }

  results_list[[sp_code]] <- route_ch
  cat("  Done:", nrow(route_ch), "routes\n")
}

scenario_manifest <- target_spp %>%
  select(Common.Name, Code, bbs_english, model_scenario, min_max_route_years,
         iter_warmup, iter_sampling, chains, adapt_delta, max_treedepth,
         any_of(c("assignment_reason", "diagnostic_flags")))
write.csv(scenario_manifest,
          file.path(fit_dir, paste0(land_cover, "_scenario_manifest.csv")),
          row.names = FALSE)

cat("\n=== Summary ===\n")
cat("Species processed this run:", length(results_list), "\n")
cat("Scenario fits in:", fit_dir, "\n")
cat("Scenario route CSVs in:", route_out_dir, "\n")
cat("Scenario manifest:", file.path(fit_dir, paste0(land_cover, "_scenario_manifest.csv")), "\n")
