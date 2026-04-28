# Adapted from AdamCSmithCWS/Route-level_BBS_trends/1alt_Species_data_prep_bbsBayes2.R
# & AdamCSmithCWS/Jefferys_etal/Species_data_prep_mean_habitat_and_slope.R
# null model adapted from AdamCSmithCWS/Jefferys_etal/Fitting_new_iCAR_slope_model.R
#
# Generic version: set `target_name` below to switch between covariates
# (e.g. "developed", "grassland"). Requires the covariate CSV at
#   data/<target_name>_1km/<target_name>.csv
# with the same schema as the developed/grassland files: a column named
# `<target_name>` plus StateNum, Route, year (etc.).

library(bbsBayes2)
library(tidyverse)
library(cmdstanr)
library(patchwork)
library(sf)
library(here)

here::i_am("code/1_Species_data_prep_covariate_mean_and_slope.R")
source("functions/neighbours_define_voronoi.R") ## function to define neighbourhood relationships for spatial model comparison


# User settings -----------------------------------------------------------
strat   <- "bbs_usgs"
model   <- "slope"
species <- "Baird's Sparrow"

firstYear <- 2010
lastYear  <- 2024
year_range <- firstYear:lastYear

# Covariate to use as the route-level habitat predictor. Must match the
# sub-folder name under data/ (with "_1km" suffix), the CSV filename, and
# the covariate column name inside that CSV. E.g. "developed" or "grassland".
target_name <- "grassland"

# Derived names used throughout
cov_dir      <- paste0(target_name, "_1km")                       # e.g. "grassland_1km"
cov_csv      <- here::here("data", cov_dir, paste0(target_name, ".csv"))
mean_col     <- paste0("mean_",  target_name)                     # e.g. "mean_grassland"
slope_col    <- paste0("slope_", target_name)                     # e.g. "slope_grassland"


# BBS data ----------------------------------------------------------------
data_pkg <- bbsBayes2::stratify(by = strat, species = species,
                                use_map = FALSE) %>%
  bbsBayes2::prepare_data(min_year = firstYear,
                          max_year = lastYear,
                          min_n_routes = 1,
                          min_max_route_years = 1)

raw_data <- data_pkg[["raw_data"]] # package now retains the lat long information for each route

# Include only continental US
raw_data <- raw_data %>%
  filter(country_num == 840) %>%
  filter(state_num != 3)

# strata map of one of the bbsBayes base maps
# helps group and set boundaries for the route-level neighbours,
## NOT directly used in the model
strata_map <- load_map(strat)

# create list of routes and locations to ID routes that are not inside of original
# strata (some off-shore islands)
route_map1 <- raw_data %>%
  select(route, strata_name, latitude, longitude) %>%
  distinct()

route_map1 <- st_as_sf(route_map1, coords = c("longitude", "latitude"))
st_crs(route_map1) <- 4326 # BBS database indicates that coordinates are stored in WGS 84
route_map1 <- st_transform(route_map1, crs = st_crs(strata_map))

# drops routes geographically outside of the strata (some offshore islands)
# and adds the strat indicator variable to link to model output
strata_map_buf <- strata_map %>%
  filter(strata_name %in% route_map1$strata_name) %>%
  summarise() %>%
  st_buffer(., 10000) # drops any routes with start-points > 10 km outside of strata boundaries
realized_routes <- route_map1 %>%
  st_join(., strata_map_buf,
          join = st_within,
          left = FALSE)

# reorganizes data after routes were dropped outside of strata
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


# Load covariate ----------------------------------------------------------
# Load BEFORE the spatial/neighbourhood setup so that routes with no covariate
# data are dropped before routeF indices and the adjacency graph are built.
# Otherwise stan_data[["nroutes"]], the adjacency, and the habitat vectors
# end up with mismatched dimensions.
cov_raw <- read.csv(cov_csv) %>%
  mutate(rt.uni = paste(StateNum, Route, sep = "-"))

# record routes missing covariate information (optional, commented out)
# missing_cov <- cov_raw %>% filter(is.na(.data[[target_name]]))
# if(nrow(missing_cov) > 0){
#   if(!dir.exists(here::here("data","missing_covariates"))) dir.create(here::here("data","missing_covariates"), recursive = TRUE)
#   write.csv(missing_cov,
#             here::here("data","missing_covariates",
#                        paste0("routes_missing_",target_name,"_",firstYear,"_",lastYear,".csv")),
#             row.names = FALSE)
# }

# Tidy covariate table with a single generic column `cov_value` so the rest of
# the script doesn't need to know the original column name.
cov_full <- cov_raw %>%
  filter(!is.na(.data[[target_name]])) %>%
  rename(route = rt.uni) %>%
  transmute(route,
            year,
            cov_value = .data[[target_name]]) %>%
  filter(year %in% year_range)

# drop count data (and routes) with no covariate information
new_data <- new_data %>%
  inner_join(cov_full,
             by = c("route", "r_year" = "year"))

strata_list <- data.frame(strata_name = unique(new_data$strat_name),
                          strat       = unique(new_data$strat))

realized_strata_map <- strata_map %>%
  filter(strata_name %in% strata_list$strata_name)


# Spatial neighbours set up -----------------------------------------------
new_data$routeF <- as.integer(factor((new_data$route))) # main route-level integer index

# create a data frame of each unique route in the species-specific dataset
route_map <- unique(data.frame(route     = new_data$route,
                               routeF    = new_data$routeF,
                               strat     = new_data$strat_name,
                               latitude  = new_data$latitude,
                               longitude = new_data$longitude))

# reconcile duplicate spatial locations
# adhoc way of separating different routes with the same starting coordinates
# this shifts the starting coordinates of the duplicates by ~1.5km to the North East
# ensures that the duplicates have a unique spatial location, but remain very close to
# their original location and retain reasonable neighbourhood relationships
dups <- which(duplicated(route_map[, c("latitude", "longitude")]))
while (length(dups) > 0) {
  route_map[dups, "latitude"]  <- route_map[dups, "latitude"]  + 0.01 # ~1km
  route_map[dups, "longitude"] <- route_map[dups, "longitude"] + 0.01 # ~1km
  dups <- which(duplicated(route_map[, c("latitude", "longitude")]))
}
if (length(which(duplicated(route_map[, c("latitude", "longitude")]))) > 0) {
  stop("At least one duplicate route remains")
}

# create spatial object from route_map dataframe
route_map <- st_as_sf(route_map, coords = c("longitude", "latitude"))
st_crs(route_map) <- 4326 # WGS 84
route_map <- st_transform(route_map, crs = st_crs(strata_map))

car_stan_dat <- neighbours_define_voronoi(real_point_map = route_map,
                                          species        = species,
                                          strat_indicator = "routeF",
                                          strata_map     = realized_strata_map,
                                          concavity      = 1) # concavity arg from concaveman()

print(car_stan_dat$map)

# save the map of connections to a PDF for inspection
species_f <- gsub(gsub(species, pattern = " ", replacement = "_", fixed = TRUE),
                  pattern = "'", replacement = "", fixed = TRUE)
if (!dir.exists(here::here("data", "maps"))) dir.create(here::here("data", "maps"), recursive = TRUE)
pdf(here::here("data", "maps",
               paste0("route_map_", firstYear, "-", lastYear, "_", species_f, "_", target_name, ".pdf")))
print(car_stan_dat$map)
dev.off()


# Route-level mean and slope of the covariate ------------------------------
route_df <- route_map %>% sf::st_drop_geometry()

sl <- function(y, x) {
  t <- lm(y ~ x)
  coefficients(t)[["x"]]
}

# restrict covariate summaries to the routes that survived the BBS/adjacency
# filtering so slope_full lines up 1:1 with routeF 1..nroutes
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

# sanity check: one row per route, ordered by routeF
stopifnot(nrow(slope_full) == max(new_data$routeF))
stopifnot(identical(slope_full$routeF, seq_len(nrow(slope_full))))


# Build the data list required for Stan -----------------------------------
stan_data <- list()
stan_data[["count"]]     <- new_data$count
stan_data[["ncounts"]]   <- length(new_data$count)
stan_data[["strat"]]     <- new_data$strat
stan_data[["route"]]     <- new_data$routeF
stan_data[["year"]]      <- new_data$year
stan_data[["firstyr"]]   <- new_data$firstyr
stan_data[["fixedyear"]] <- floor(stats::median(year_range)) # mid-year of time series

stan_data[["nyears"]]     <- max(new_data$year)
stan_data[["observer"]]   <- as.integer(factor(new_data$ObsN))
stan_data[["nobservers"]] <- max(stan_data$observer)

stan_data[["N_edges"]] <- car_stan_dat$N_edges
stan_data[["node1"]]   <- car_stan_dat$node1
stan_data[["node2"]]   <- car_stan_dat$node2
stan_data[["nroutes"]] <- max(stan_data$route)

# Predictor transformation follows Smith et al. 2024 (ACE 19(2):2) Appendix 2 /
# Jefferys et al. Rufous Hummingbird example:
#   - route_habitat: mean-center the route-level mean habitat value only.
#     Units stay as "proportion of <target_name> above the range-wide average",
#     matching the prior rho_ALPHA_hab ~ N(0,1).
#   - route_habitat_slope: mean-center AND multiply by 100, so the slope is in
#     percentage-points per year above the range-wide average, matching the
#     priors rho_BETA_hab ~ N(0,0.1) and sdrho_beta_hab ~ N(0,0.1).
# Avoid scale(): SD-scaling makes prior strength depend on species/time-window.
mean_vec  <- slope_full[[mean_col]]
slope_vec <- slope_full[[slope_col]]

stan_data[["route_habitat"]]       <- as.numeric(mean_vec  - mean(mean_vec,  na.rm = TRUE))
stan_data[["route_habitat_slope"]] <- 100 * (slope_vec - mean(slope_vec, na.rm = TRUE))

if (car_stan_dat$N != stan_data[["nroutes"]]) stop("Some routes are missing from adjacency matrix")

dist_matrix_km <- dist_matrix(route_map, strat_indicator = "routeF")


# Save prepped data so model fitting can be re-run without re-preparing ----
if (!dir.exists(here::here("data", "stan_data"))) dir.create(here::here("data", "stan_data"), recursive = TRUE)
sp_data_file <- here::here("data", "stan_data",
                           paste0(species_f, "_", target_name, "_",
                                  firstYear, "_", lastYear, "_stan_data.RData"))
save(list = c("stan_data",
              "new_data",
              "route_map",
              "realized_strata_map",
              "car_stan_dat",
              "dist_matrix_km",
              "cov_full",
              "slope_full",
              "target_name"),
     file = sp_data_file)


# Fit ---------------------------------------------------------------------
## for both time-periods, there is relatively strong spatial autocorrelation
## in both the habitat suitability and the mean abundance of the species.
## Since the spatial component of habitat suitability could reasonably be a
## cause of the spatial dependency in abundance, we estimated the residual
## component of the intercept with a non-spatial (simple random effect).
## Setting `spatial_intercept <- TRUE` switches the residual intercept back to
## the spatial (iCAR) term.
spatial_intercept <- TRUE
stan_data[["fit_spatial"]] <- ifelse(spatial_intercept, 1, 0)

mod.file <- "models/slope_habitat_route_NB.stan"
slope_model <- cmdstan_model(mod.file, stanc_options = list("Oexperimental"))

stanfit <- slope_model$sample(
  data             = stan_data,
  refresh          = 400,
  iter_sampling    = 2000,
  iter_warmup      = 2000,
  max_treedepth    = 15,
  parallel_chains  = 4,
  save_cmdstan_config = TRUE) # needed for $summary() in cmdstanr >= 0.8.0

summ <- stanfit$summary()


# Save fit ----------------------------------------------------------------
output_dir <- "output"
if (!dir.exists(output_dir)) dir.create(output_dir)

out_base <- paste0(species_f, "_", target_name, "_", firstYear, "_", lastYear)

print(paste(species, stanfit$time()[["total"]]))

saveRDS(stanfit, file.path(output_dir, paste0(out_base, "_stanfit.rds")))
saveRDS(summ,    file.path(output_dir, paste0(out_base, "_summ_fit.rds")))

summ %>% arrange(-rhat)
summ %>% filter(variable %in% c("BETA", "rho_BETA_hab"))
summ %>% filter(variable %in% c("ALPHA", "rho_ALPHA_hab"))
summ %>% filter(grepl("T", variable))
summ %>% filter(grepl("CH", variable))


# Map trends --------------------------------------------------------------
mn0 <- new_data %>%
  group_by(routeF) %>%
  summarise(mn = mean(count),
            mx = max(count),
            ny = n(),
            fy = min(year),
            ly = max(year),
            sp = max(year) - min(year))

exp_t <- function(x) (exp(x) - 1) * 100

base_strata_map <- bbsBayes2::load_map("bbs_usgs")

strata_bounds <- st_union(route_map)
bb <- st_bbox(strata_bounds)
xlms <- as.numeric(c(bb$xmin, bb$xmax))
ylms <- as.numeric(c(bb$ymin, bb$ymax))

betas1 <- summ %>%
  filter(grepl("beta[", variable, fixed = TRUE)) %>%
  mutate(across(2:7, ~exp_t(.x)),
         routeF    = as.integer(str_extract(variable, "[[:digit:]]{1,}")),
         parameter = "Full with Habitat-Change") %>%
  select(routeF, mean, sd, parameter) %>%
  rename(trend = mean, trend_se = sd)

alpha1 <- summ %>%
  filter(grepl("alpha[", variable, fixed = TRUE)) %>%
  mutate(across(2:7, ~exp(.x)),
         routeF    = as.integer(str_extract(variable, "[[:digit:]]{1,}")),
         parameter = "Full with Habitat") %>%
  select(routeF, median, sd) %>%
  rename(abundance = median, abundance_se = sd)

alpha2 <- summ %>%
  filter(grepl("alpha_resid[", variable, fixed = TRUE)) %>%
  mutate(across(2:7, ~exp(.x)),
         routeF    = as.integer(str_extract(variable, "[[:digit:]]{1,}")),
         parameter = "Residual") %>%
  select(routeF, median, sd) %>%
  rename(abundance = median, abundance_se = sd)

betas1 <- betas1 %>% inner_join(alpha1)

betas2 <- summ %>%
  filter(grepl("beta_resid[", variable, fixed = TRUE)) %>%
  mutate(across(2:7, ~exp_t(.x)),
         routeF    = as.integer(str_extract(variable, "[[:digit:]]{1,}")),
         parameter = "Residual") %>%
  select(routeF, mean, sd, parameter) %>%
  rename(trend = mean, trend_se = sd)
betas2 <- betas2 %>%
  inner_join(alpha2, by = "routeF") %>%
  inner_join(mn0,    by = "routeF")

betas <- bind_rows(betas1, betas2)

plot_map <- route_map %>%
  left_join(betas, by = "routeF", multiple = "all")

breaks <- c(-7, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 7)
lgnd_head   <- "Mean Trend\n"
trend_title <- "Mean Trend"
labls <- c(paste0("< ", breaks[1]),
           paste0(breaks[-length(breaks)], ":", breaks[-1]),
           paste0("> ", breaks[length(breaks)]))
labls <- paste0(labls, " %/year")
plot_map$Tplot <- cut(plot_map$trend, breaks = c(-Inf, breaks, Inf), labels = labls)

map_palette <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf",
                 "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")
names(map_palette) <- labls

scbar <- scalebar(plot_map,
                  dist      = 300,
                  dist_unit = "km",
                  transform = FALSE,
                  facet.var = "parameter",
                  facet.lev = "Full with Habitat-Change",
                  location  = "bottomleft",
                  st.size   = 2.5,
                  box.fill  = c(gray(0.5), "white"),
                  box.color = gray(0.3),
                  st.color  = gray(0.5))

map <- ggplot() +
  geom_sf(data = base_strata_map, fill = NA, colour = grey(0.75)) +
  geom_sf(data = plot_map, aes(colour = Tplot, size = abundance)) +
  scale_size_continuous(range = c(0.05, 2), name = "Mean Count") +
  scale_colour_manual(values = map_palette, aesthetics = c("colour"),
                      guide  = guide_legend(reverse = TRUE),
                      name   = lgnd_head) +
  coord_sf(xlim = xlms, ylim = ylms) +
  guides(size = "none") +
  xlab("") + ylab("") +
  labs(title = "Trend") +
  theme_bw() +
  facet_wrap(vars(parameter))

map <- map + scbar
map

map_abund <- ggplot() +
  geom_sf(data = base_strata_map, fill = NA, colour = grey(0.75)) +
  geom_sf(data = plot_map, aes(colour = abundance)) +
  scale_colour_viridis_c(begin = 0.1, end = 0.9,
                         guide = guide_legend(reverse = TRUE),
                         name  = "Relative Abundance") +
  coord_sf(xlim = xlms, ylim = ylms) +
  theme_bw() +
  xlab("") + ylab("") +
  labs(title = "Relative Abundance") +
  facet_wrap(vars(parameter))

map_se <- ggplot() +
  geom_sf(data = base_strata_map, fill = NA, colour = grey(0.75)) +
  geom_sf(data = plot_map, aes(colour = trend_se, size = abundance_se)) +
  scale_size_continuous(range = c(0.05, 2),
                        name  = "SE of Mean Count",
                        trans = "reverse") +
  scale_colour_viridis_c(aesthetics = c("colour"),
                         guide = guide_legend(reverse = TRUE),
                         name  = "SE of Trend") +
  coord_sf(xlim = xlms, ylim = ylms) +
  theme_bw() +
  xlab("") + ylab("") +
  guides(size = "none") +
  labs(title = "") +
  facet_wrap(vars(parameter))
