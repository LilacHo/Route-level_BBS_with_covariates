# Calculate route-level CH_no_habitat (habitat covariate model) for all
# species in a target group (e.g., grasslands, aridlands). Outputs one CSV
# per species with:
#   species, route, routeF, latitude, longitude, CH, CH_no_habitat, CH_dif
#
# CH_no_habitat_route[s] = 100 * (exp(beta_resid[s] * dt) - 1)
#   where dt = lastYear - firstYear, and beta_resid[s] is the posterior mean
#   of the residual slope (trend net of the habitat-change component for the
#   target group).
#
# If pre-fitted output files (summ_fit.rds + stan_data.RData) do not exist for
# a species, the full data-prep and model-fitting pipeline is run inline to
# produce them.
#
# To run for a different group, change `target_name` in the Settings block.

# install.packages("concaveman")
library(bbsBayes2)
# fetch_bbs_data()
library(tidyverse)
library(cmdstanr)
library(posterior)
library(sf)
library(here)

here::i_am("code/1_CH_no_habitat_routes.R")
source("functions/neighbours_define_voronoi.R")

# Settings ----------------------------------------------------------------
target_name <- "grasslands"   # target group: must match a value in the
                              # 'Group' column of spp_names_codes_group_aou.csv
                              # and the covariate file data/<target_name>.csv
firstYear   <- 2010
lastYear    <- 2024
year_range  <- firstYear:lastYear
dt          <- lastYear - firstYear   # 14

strat <- "bbs_usgs"
spatial_intercept <- TRUE

# Derived covariate names
cov_csv  <- here::here("data",  paste0(target_name, ".csv"))
mean_col <- paste0("mean_", target_name)
slope_col <- paste0("slope_", target_name)

# Target group species list -----------------------------------------------
spp_df <- read.csv(here::here("data", "spp_names_codes_group_aou.csv"))

target_spp <- spp_df %>%
  filter(Group == target_name) %>%
  distinct(Common.Name, Code, .keep_all = TRUE)

cat(target_name, "species (n =", nrow(target_spp), "):\n")
print(target_spp %>% select(Common.Name, Code, bbs_english))

# Helper: convert species name to file-safe format (matches data-prep script)
species_to_f <- function(sp) {
  gsub("'", "", gsub(" ", "_", sp, fixed = TRUE), fixed = TRUE)
}

# Compile the Stan model once (reused across species) ---------------------
mod.file <- "models/slope_habitat_route_NB.stan"
slope_model <- cmdstan_model(mod.file, force_recompile = TRUE)

# Helper: slope of a linear regression
sl <- function(y, x) {
  t <- lm(y ~ x)
  coefficients(t)[["x"]]
}

# CH_no_habitat calculation helper
pct_cum <- function(b, dt) 100 * (exp(b * dt) - 1)

# ==========================================================================
# Function: run the full data-prep + fit pipeline for one species
# ==========================================================================
fit_species <- function(species, species_bbs, species_f, target_name, strat, firstYear, lastYear,
                        year_range, cov_csv, mean_col, slope_col,
                        spatial_intercept, slope_model) {

  cat("    Running full data-prep + model fit for:", species, "\n")

  # BBS data ---------------------------------------------------------------
  data_pkg <- bbsBayes2::stratify(by = strat, species = species_bbs,
                                  use_map = FALSE) %>%
    bbsBayes2::prepare_data(min_year = firstYear,
                            max_year = lastYear,
                            min_n_routes = 1,
                            min_max_route_years = 1)

  raw_data <- data_pkg[["raw_data"]]

  # Continental US only

  raw_data <- raw_data %>%
    filter(country_num == 840) %>%
    filter(state_num != 3)

  strata_map <- load_map(strat)

  # Route filtering --------------------------------------------------------
  route_map1 <- raw_data %>%
    select(route, strata_name, latitude, longitude) %>%
    distinct()

  route_map1 <- st_as_sf(route_map1, coords = c("longitude", "latitude"))
  st_crs(route_map1) <- 4326
  route_map1 <- st_transform(route_map1, crs = st_crs(strata_map))

  strata_map_buf <- strata_map %>%
    filter(strata_name %in% route_map1$strata_name) %>%
    summarise() %>%
    st_buffer(., 10000)
  realized_routes <- route_map1 %>%
    st_join(., strata_map_buf, join = st_within, left = FALSE)

  new_data <- data.frame(strat_name = raw_data$strata_name,
                         strat      = raw_data$strata,
                         route      = raw_data$route,
                         latitude   = raw_data$latitude,
                         longitude  = raw_data$longitude,
                         count      = raw_data$count,
                         year       = raw_data$year_num,
                         firstyr    = raw_data$first_year,
                         ObsN       = raw_data$observer,
                         r_year     = raw_data$year) %>%
    filter(route %in% realized_routes$route)

  # Load covariate ---------------------------------------------------------
  cov_raw <- read.csv(cov_csv) %>%
    mutate(rt.uni = paste(StateNum, Route, sep = "-"))

  cov_full <- cov_raw %>%
    filter(!is.na(.data[[target_name]])) %>%
    rename(route = rt.uni) %>%
    transmute(route,
              year,
              cov_value = .data[[target_name]]) %>%
    filter(year %in% year_range)

  # Drop routes with no covariate
  new_data <- new_data %>%
    inner_join(cov_full, by = c("route", "r_year" = "year"))

  strata_list <- data.frame(strata_name = unique(new_data$strat_name),
                            strat       = unique(new_data$strat))
  realized_strata_map <- strata_map %>%
    filter(strata_name %in% strata_list$strata_name)

  # Spatial neighbours -----------------------------------------------------
  new_data$routeF <- as.integer(factor((new_data$route)))

  route_map <- unique(data.frame(route     = new_data$route,
                                 routeF    = new_data$routeF,
                                 strat     = new_data$strat_name,
                                 latitude  = new_data$latitude,
                                 longitude = new_data$longitude))

  # Fix duplicate coordinates
  dups <- which(duplicated(route_map[, c("latitude", "longitude")]))
  while (length(dups) > 0) {
    route_map[dups, "latitude"]  <- route_map[dups, "latitude"]  + 0.01
    route_map[dups, "longitude"] <- route_map[dups, "longitude"] + 0.01
    dups <- which(duplicated(route_map[, c("latitude", "longitude")]))
  }
  if (length(which(duplicated(route_map[, c("latitude", "longitude")]))) > 0) {
    stop("At least one duplicate route remains")
  }

  route_map <- st_as_sf(route_map, coords = c("longitude", "latitude"))
  st_crs(route_map) <- 4326
  route_map <- st_transform(route_map, crs = st_crs(strata_map))

  car_stan_dat <- neighbours_define_voronoi(real_point_map = route_map,
                                            species        = species,
                                            strat_indicator = "routeF",
                                            strata_map     = realized_strata_map,
                                            concavity      = 1)

  # Save map PDF
  if (!dir.exists(here::here("data", "maps"))) dir.create(here::here("data", "maps"), recursive = TRUE)
  pdf(here::here("data", "maps",
                 paste0("route_map_", firstYear, "-", lastYear, "_", species_f, "_", target_name, ".pdf")))
  print(car_stan_dat$map)
  dev.off()

  # Route-level mean and slope of covariate --------------------------------
  route_df <- route_map %>% sf::st_drop_geometry()
  cov_full <- cov_full %>% filter(route %in% route_df$route)

  cov_mean_t <- cov_full %>%
    group_by(route) %>%
    summarise(!!mean_col := mean(cov_value, na.rm = TRUE), .groups = "drop")

  cov_slope_t <- cov_full %>%
    arrange(year, route) %>%
    group_by(route) %>%
    summarise(!!slope_col := sl(cov_value, year), .groups = "drop")

  slope_full <- inner_join(cov_mean_t, cov_slope_t, by = "route") %>%
    mutate(startYear = firstYear) %>%
    left_join(route_df, by = "route") %>%
    arrange(routeF) %>%
    filter(!is.na(routeF))

  stopifnot(nrow(slope_full) == max(new_data$routeF))
  stopifnot(identical(slope_full$routeF, seq_len(nrow(slope_full))))

  # Build Stan data --------------------------------------------------------
  stan_data <- list()
  stan_data[["count"]]     <- new_data$count
  stan_data[["ncounts"]]   <- length(new_data$count)
  stan_data[["strat"]]     <- new_data$strat
  stan_data[["route"]]     <- new_data$routeF
  stan_data[["year"]]      <- new_data$year
  stan_data[["firstyr"]]   <- new_data$firstyr
  stan_data[["fixedyear"]] <- floor(stats::median(year_range))

  stan_data[["nyears"]]     <- max(new_data$year)
  stan_data[["observer"]]   <- as.integer(factor(new_data$ObsN))
  stan_data[["nobservers"]] <- max(stan_data$observer)

  stan_data[["N_edges"]] <- car_stan_dat$N_edges
  stan_data[["node1"]]   <- car_stan_dat$node1
  stan_data[["node2"]]   <- car_stan_dat$node2
  stan_data[["nroutes"]] <- max(stan_data$route)

  mean_vec  <- slope_full[[mean_col]]
  slope_vec <- slope_full[[slope_col]]

  stan_data[["route_habitat"]]       <- as.numeric(mean_vec - mean(mean_vec, na.rm = TRUE))
  stan_data[["route_habitat_slope"]] <- 100 * (slope_vec - mean(slope_vec, na.rm = TRUE))

  if (car_stan_dat$N != stan_data[["nroutes"]]) stop("Some routes are missing from adjacency matrix")

  dist_matrix_km <- dist_matrix(route_map, strat_indicator = "routeF")

  stan_data[["fit_spatial"]] <- ifelse(spatial_intercept, 1, 0)

  # Save stan_data ---------------------------------------------------------
  if (!dir.exists(here::here("data", "stan_data"))) dir.create(here::here("data", "stan_data"), recursive = TRUE)
  sp_data_file <- here::here("data", "stan_data",
                             paste0(species_f, "_", target_name, "_",
                                    firstYear, "_", lastYear, "_stan_data.RData"))
  save(list = c("stan_data", "new_data", "route_map", "realized_strata_map",
                "car_stan_dat", "dist_matrix_km", "cov_full", "slope_full", "target_name"),
       file = sp_data_file)

  # Fit model --------------------------------------------------------------
  stanfit <- slope_model$sample(
    data             = stan_data,
    refresh          = 400,
    chains           = 4,
    parallel_chains  = 4,
    iter_sampling    = 2000,
    iter_warmup      = 2000,
    max_treedepth    = 15,
    show_exceptions  = TRUE,
    save_cmdstan_config = TRUE)

  summ <- stanfit$summary()

  # Save fit ---------------------------------------------------------------
  output_dir <- here::here("output")
  if (!dir.exists(output_dir)) dir.create(output_dir)

  out_base <- paste0(species_f, "_", target_name, "_", firstYear, "_", lastYear)
  saveRDS(stanfit, file.path(output_dir, paste0(out_base, "_stanfit.rds")))
  saveRDS(summ,    file.path(output_dir, paste0(out_base, "_summ_fit.rds")))

  cat("    Fit complete for:", species, "(", stanfit$time()[["total"]], "s )\n")

  # Return summary and route_map
  list(summ = summ, route_map = route_map)
}

# ==========================================================================
# Main loop: compute CH_no_habitat per route for each species in the group
# ==========================================================================
results_list <- list()

for (i in seq_len(nrow(target_spp))) {
  sp      <- target_spp$Common.Name[i]
  sp_f    <- species_to_f(sp)
  sp_code <- target_spp$Code[i]
  sp_bbs  <- target_spp$bbs_english[i]   # name used to look up species in BBS

  cat("\n[", i, "/", nrow(target_spp), "]", sp, "\n")

  # Skip if the per-species output CSV already exists
  sp_out_dir <- here::here("output", "species_routes")
  sp_csv <- file.path(sp_out_dir, paste0(sp_f, "_route_CH.csv"))
  if (file.exists(sp_csv)) {
    cat("  Skipping (already exists):", basename(sp_csv), "\n")
    next
  }

  # Paths to pre-fitted output

  summ_file <- here::here("output",
                          paste0(sp_f, "_", target_name, "_",
                                 firstYear, "_", lastYear, "_summ_fit.rds"))
  stan_data_file <- here::here("data", "stan_data",
                               paste0(sp_f, "_", target_name, "_",
                                      firstYear, "_", lastYear, "_stan_data.RData"))

  # If both files exist, load them; otherwise run full pipeline
  if (file.exists(summ_file) && file.exists(stan_data_file)) {
    cat("  Loading existing fit\n")
    summ <- readRDS(summ_file)
    load(stan_data_file)  # loads route_map (sf) and other objects
  } else {
    cat("  No existing fit found — running data prep + model fit\n")
    fit_result <- tryCatch(
      fit_species(species = sp,
                  species_bbs = sp_bbs,
                  species_f = sp_f,
                  target_name = target_name,
                  strat = strat,
                  firstYear = firstYear,
                  lastYear = lastYear,
                  year_range = year_range,
                  cov_csv = cov_csv,
                  mean_col = mean_col,
                  slope_col = slope_col,
                  spatial_intercept = spatial_intercept,
                  slope_model = slope_model),
      error = function(e) {
        message("  [ERROR] Failed for ", sp, ": ", conditionMessage(e))
        # Try to get CmdStan output files for diagnosis
        out_files <- list.files(tempdir(), pattern = paste0(".*", species_to_f(sp), ".*"),
                                full.names = TRUE, recursive = TRUE)
        if (length(out_files) > 0) {
          message("  Temp files found: ", paste(out_files, collapse = "\n    "))
        }
        return(NULL)
      }
    )
    if (is.null(fit_result)) next
    summ      <- fit_result$summ
    route_map <- fit_result$route_map
  }

  # Extract beta (full slope) and beta_resid (residual slope) per route
  beta_df <- summ %>%
    filter(str_detect(variable, "^beta\\[")) %>%
    transmute(routeF = as.integer(str_extract(variable, "\\d+")),
              beta = mean)

  beta_resid_df <- summ %>%
    filter(str_detect(variable, "^beta_resid\\[")) %>%
    transmute(routeF = as.integer(str_extract(variable, "\\d+")),
              beta_resid = mean)

  # Get lat/lon from route_map (sf object)
  route_info <- route_map %>%
    st_transform(4326) %>%
    mutate(longitude = st_coordinates(.)[, 1],
           latitude  = st_coordinates(.)[, 2]) %>%
    st_drop_geometry() %>%
    select(route, routeF, latitude, longitude)

  # Compute route-level CH, CH_no_habitat, CH_dif
  route_ch <- beta_df %>%
    left_join(beta_resid_df, by = "routeF") %>%
    mutate(CH            = pct_cum(beta, dt),
           CH_no_habitat = pct_cum(beta_resid, dt),
           CH_dif        = CH - CH_no_habitat) %>%
    left_join(route_info, by = "routeF") %>%
    mutate(species      = sp,
           species_code = sp_code) %>%
    select(species, species_code, route, routeF, latitude, longitude,
           CH, CH_no_habitat, CH_dif)

  results_list[[sp]] <- route_ch

  # Save per-species CSV
  if (!dir.exists(sp_out_dir)) dir.create(sp_out_dir, recursive = TRUE)
  write.csv(route_ch,
            file.path(sp_out_dir, paste0(sp_f, "_route_CH.csv")),
            row.names = FALSE)

  cat("  Done:", nrow(route_ch), "routes\n")
}

cat("\n=== Summary ===\n")
cat("Species processed this run:", length(results_list), "\n")
cat("Per-species CSVs in:", here::here("output", "species_routes"), "\n")
