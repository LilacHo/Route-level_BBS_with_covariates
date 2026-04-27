# Adapted from AdamCSmithCWS/Route-level_BBS_trends/1alt_Species_data_prep_bbsBayes2.R 
# & AdamCSmithCWS/Jefferys_etal/Species_data_prep_mean_habitat_and_slope.R
# null model adapted from AdamCSmithCWS/Jefferys_etal/Fitting_new_iCAR_slope_model.R

library(bbsBayes2)
library(tidyverse)
library(cmdstanr)
library(patchwork)
library(sf)
library(here)

here::i_am("code/1 Species_data_prep_developed_mean_and_slope.R")
source("functions/neighbours_define_voronoi.R") ## function to define neighbourhood relationships for spatial model comparison


strat = "bbs_usgs"
model = "slope"
species = "Baird's Sparrow"

firstYear <- 2010
lastYear <- 2024
year_range <- c(2010:2024)


data_pkg <- bbsBayes2::stratify(by = strat,species = species,
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


#create list of routes and locations to ID routes that are not inside of original strata (some off-shore islands)
route_map1 <- raw_data %>% 
  select(route,strata_name,latitude,longitude) %>% 
  distinct()
# unique(data.frame(route = raw_data$route,
#                             strat = raw_data$strata_name,
#                             latitude = raw_data$latitude,
#                             longitude = raw_data$longitude))

#create spatial object of above
route_map1 <- st_as_sf(route_map1,coords = c("longitude","latitude"))
st_crs(route_map1) <- 4326 #BBS database indicates that coordinates are now stored in WGS 84

#reconcile the projections of routes and base bbs strata
route_map1 = st_transform(route_map1,crs = st_crs(strata_map))

#drops the routes geographically outside of the strata (some offshore islands) 
# and adds the strat indicator variable to link to model output
strata_map_buf <- strata_map %>% 
  filter(strata_name %in% route_map1$strata_name) %>% 
  summarise() %>% 
  st_buffer(.,10000) #drops any routes with start-points > 10 km outside of strata boundaries
realized_routes <- route_map1 %>% 
  st_join(.,strata_map_buf,
          join = st_within,
          left = FALSE) 



# reorganizes data after routes were dropped outside of strata
new_data <- data.frame(strat_name = raw_data$strata_name,
                       strat = raw_data$strata,
                       route = raw_data$route,
                       latitude = raw_data$latitude,
                       longitude = raw_data$longitude,
                       count = raw_data$count,
                       year = raw_data$year_num,
                       firstyr = raw_data$first_year,
                       ObsN = raw_data$observer,
                       r_year = raw_data$year) %>% 
  filter(route %in% realized_routes$route)

strata_list <- data.frame(strata_name = unique(new_data$strat_name),
                          strat = unique(new_data$strat))


realized_strata_map <- strata_map %>%
  filter(strata_name %in% strata_list$strata_name)

# Spatial neighbours set up --------------------

new_data$routeF <- as.integer(factor((new_data$route))) #main route-level integer index

#create a data frame of each unique route in the species-specific dataset
route_map = unique(data.frame(route = new_data$route,
                              routeF = new_data$routeF,
                              strat = new_data$strat_name,
                              latitude = new_data$latitude,
                              longitude = new_data$longitude))


# reconcile duplicate spatial locations -----------------------------------
# adhoc way of separating different routes with the same starting coordinates
# this shifts the starting coordinates of the duplicates by ~1.5km to the North East 
# ensures that the duplicates have a unique spatial location, but remain very close to
# their original location and retain reasonable neighbourhood relationships
# these duplicates happen when a "new" route is established (i.e., a route is re-named) 
# because some large proportion
# of the end of a route is changed, but the start-point remains the same
dups = which(duplicated(route_map[,c("latitude","longitude")]))
while(length(dups) > 0){
  route_map[dups,"latitude"] <- route_map[dups,"latitude"]+0.01 #=0.01 decimal degrees ~ 1km
  route_map[dups,"longitude"] <- route_map[dups,"longitude"]+0.01 #=0.01 decimal degrees ~ 1km
  dups = which(duplicated(route_map[,c("latitude","longitude")]))
  
}
dups = which(duplicated(route_map[,c("latitude","longitude")])) 
if(length(dups) > 0){stop(paste(spec,"At least one duplicate route remains"))}



#create spatial object from route_map dataframe
route_map = st_as_sf(route_map,coords = c("longitude","latitude"))
st_crs(route_map) <- 4326 #BBS database indicates that coordinates are now stored in WGS 84


#reproject teh routes spatial object ot match the strata-map in equal area projection
route_map = st_transform(route_map,crs = st_crs(strata_map))

car_stan_dat <- neighbours_define_voronoi(real_point_map = route_map,
                                          species = species,
                                          strat_indicator = "routeF",
                                          strata_map = realized_strata_map,
                                          concavity = 1)#concavity argument from concaveman()

print(car_stan_dat$map)

# stan_data <- list()
# stan_data[["count"]] <- new_data$count
# stan_data[["ncounts"]] <- length(new_data$count)
# stan_data[["strat"]] <- new_data$strat
# stan_data[["route"]] <- new_data$routeF
# stan_data[["year"]] <- new_data$year
# stan_data[["firstyr"]] <- new_data$firstyr
# stan_data[["fixedyear"]] <- floor(mean(c(1:max(new_data$year))))
# 
# 
# stan_data[["nyears"]] <- max(new_data$year)
# stan_data[["observer"]] <- as.integer(factor((new_data$ObsN)))
# stan_data[["nobservers"]] <- max(stan_data$observer)
# 
# 
# 
# stan_data[["N_edges"]] <- car_stan_dat$N_edges
# stan_data[["node1"]] <- car_stan_dat$node1
# stan_data[["node2"]] <- car_stan_dat$node2
# stan_data[["nroutes"]] <- max(stan_data$route)


## Load covariates ####
developed <- read.csv(here::here("data", "developed_1km", "developed.csv"))
developed <- developed %>%
  mutate(rt.uni = paste(StateNum, Route, sep = "-"))

# drop the routes with no covariate information
developed <- developed %>% 
  filter(!is.na(developed))

developed_full <- developed %>% 
  rename(route = rt.uni) %>% 
  select(route, year, developed) 

# model_lm <- lm(developed ~ year, data = developed_full)
# summary(model_lm)

# developed_vis <- ggplot(data = developed,
#                   aes(x = year,y = developed,
#                       group = rt.uni))+
#   geom_line(alpha = 0.2)
# developed_vis

# Calculate slope and mean -----------------------------------------
# #create a data frame of each unique route in the species-specific dataset
# route_map2 = unique(data.frame(route = new_data$route,
#                                routeF = new_data$routeF,
#                                strat = new_data$strat_name,
#                                latitude = new_data$latitude,
#                                longitude = new_data$longitude))
# 

# route_map2 <- route_map2 %>%
#   separate(route, into = c("StateNum", "RouteNum"), sep = "-", convert = TRUE)
route_df <- route_map %>% 
  sf::st_drop_geometry(route_map)


sl <- function(y,x){
  t <- lm(y~x)
  cc <- coefficients(t)[["x"]]
  return(cc)
}

slope_full <- NULL

developed_route_mean_t <- developed_full %>% 
  #filter(year < (firstYear + 10)) %>%  #first 10 years
  group_by(route) %>% 
  summarise(mean_developed = mean(developed,na.rm = T))

developed_route_slope_t <- developed_full %>% 
  arrange(year,route) %>% 
  group_by(route) %>% 
  summarise(slope_developed = sl(developed,year))

developed_tmp <- inner_join(developed_route_mean_t,developed_route_slope_t,
                      by = "route")

slope_full <- bind_rows(slope_full,developed_tmp)

slope_full <- slope_full %>% 
  left_join(.,route_df,
            by = "route") %>% 
  arrange(routeF)%>% 
  filter(!is.na(routeF))

# Drop count data with no covariate information -----------------------------------------

new_data <- new_data %>% 
  inner_join(.,developed_full,
             by = c("route",
                    "r_year" = "year"))


### Build the data list required for Stan


stan_data <- list()
stan_data[["count"]] <- new_data$count
stan_data[["ncounts"]] <- length(new_data$count)
stan_data[["strat"]] <- new_data$strat
stan_data[["route"]] <- new_data$routeF
stan_data[["year"]] <- new_data$year
stan_data[["firstyr"]] <- new_data$firstyr
stan_data[["fixedyear"]] <- floor(stats::median(year_range)) # mid-year of time series


stan_data[["nyears"]] <- max(new_data$year)
stan_data[["observer"]] <- as.integer(factor((new_data$ObsN)))
stan_data[["nobservers"]] <- max(stan_data$observer)



stan_data[["N_edges"]] <- car_stan_dat$N_edges
stan_data[["node1"]] <- car_stan_dat$node1
stan_data[["node2"]] <- car_stan_dat$node2
stan_data[["nroutes"]] <- max(stan_data$route)
stan_data[["route_habitat"]] <- as.numeric(scale(slope_full$mean_developed))
stan_data[["route_habitat_slope"]] <- as.numeric(scale(slope_full$slope_developed))


if(car_stan_dat$N != stan_data[["nroutes"]]){stop("Some routes are missing from adjacency matrix")}

dist_matrix_km <- dist_matrix(route_map,
                              strat_indicator = "routeF")
# save(list = c("stan_data",
#               "new_data",
#               "route_map",
#               "realized_strata_map",
#               "car_stan_dat",
#               "dist_matrix_km",
#               "developed_full",
#               "slope_full"),
#      file = sp_data_file)


##Fit####
## for both time-periods, there is a relatively strong spatial autocorrelation
## in both the habitat suitability and the mean abundance of the species
## Since, the spatial component of habitat suitability could reasonably be
## considered as a cause of the spatial dependency in abundance we estimated the
## residual component of the intercept term with a non-spatial (simple random effect)
## Setting this `spatial_intercept` to TRUE will fit the model with the spatial residual term
spatial_intercept <- FALSE
# trend habitat effects are not changed, but the intercept effect is
# removes the optional spatial components for intercepts 
stan_data[["fit_spatial"]] <- ifelse(spatial_intercept,1,0)

mod.file = paste0("models/slope_habitat_route_NB.stan")


slope_model <- cmdstan_model(mod.file, stanc_options = list("Oexperimental"))

stanfit <- slope_model$sample(
  data=stan_data,
  refresh=400,
  iter_sampling=2000,
  iter_warmup=2000,
  max_treedepth = 15,
  parallel_chains = 4)

summ <- stanfit$summary()


output_dir <- "output"
 
spp1 <- "developed"
spp <- paste0("_",spp1,"_")

species_f <- gsub(gsub(species,pattern = " ",replacement = "_",fixed = T),pattern = "'",replacement = "",fixed = T)
out_base <- paste0(species_f,spp,firstYear,"_",lastYear)

print(paste(species, stanfit$time()[["total"]]))
out_base <- paste0(species_f,spp,firstYear,"_",lastYear)

saveRDS(stanfit, 
        paste0(output_dir,"/",out_base,"_stanfit.rds"))

saveRDS(summ,
        paste0(output_dir,"/",out_base,"_summ_fit.rds"))

summ %>% arrange(-rhat)
summ %>% filter(variable %in% c("BETA","rho_BETA_hab"))
summ %>% filter(variable %in% c("ALPHA","rho_ALPHA_hab"))
summ %>% filter(grepl("T",variable))
summ %>% filter(grepl("CH",variable))


mn0 <- new_data %>% 
  group_by(routeF) %>% 
  summarise(mn = mean(count),
            mx = max(count),
            ny = n(),
            fy = min(year),
            ly = max(year),
            sp = max(year)-min(year))

# route_map_2006 <- route_map 

exp_t <- function(x){
  y <- (exp(x)-1)*100
}


# plot trends -------------------------------------------------------------


base_strata_map <- bbsBayes2::load_map("bbs_usgs")


strata_bounds <- st_union(route_map) #union to provide a simple border of the realised strata
bb = st_bbox(strata_bounds)
xlms = as.numeric(c(bb$xmin,bb$xmax))
ylms = as.numeric(c(bb$ymin,bb$ymax))

betas1 <- summ %>% 
  filter(grepl("beta[",variable,fixed = TRUE)) %>% 
  mutate(across(2:7,~exp_t(.x)),
         routeF = as.integer(str_extract(variable,"[[:digit:]]{1,}")),
         parameter = "Full with Habitat-Change") %>% 
  select(routeF,mean,sd,parameter) %>% 
  rename(trend = mean,
         trend_se = sd)

alpha1 <- summ %>% 
  filter(grepl("alpha[",variable,fixed = TRUE)) %>% 
  mutate(across(2:7,~exp(.x)),
         routeF = as.integer(str_extract(variable,"[[:digit:]]{1,}")),
         parameter = "Full with Habitat") %>% 
  select(routeF,median,sd) %>% 
  rename(abundance = median,
         abundance_se = sd)

alpha2 <- summ %>% 
  filter(grepl("alpha_resid[",variable,fixed = TRUE)) %>% 
  mutate(across(2:7,~exp(.x)),
         routeF = as.integer(str_extract(variable,"[[:digit:]]{1,}")),
         parameter = "Residual") %>% 
  select(routeF,median,sd) %>% 
  rename(abundance = median,
         abundance_se = sd)

betas1 <- betas1 %>% 
  inner_join(.,alpha1)

betas2 <- summ %>% 
  filter(grepl("beta_resid[",variable,fixed = TRUE)) %>% 
  mutate(across(2:7,~exp_t(.x)),
         routeF = as.integer(str_extract(variable,"[[:digit:]]{1,}")),
         parameter = "Residual") %>% 
  select(routeF,mean,sd,parameter) %>% 
  rename(trend = mean,
         trend_se = sd)
betas2 <- betas2 %>% 
  inner_join(.,alpha2,by = "routeF") %>% 
  inner_join(.,mn0,by = "routeF")


betas <- bind_rows(betas1,betas2)

plot_map <- route_map %>% 
  left_join(.,betas,
            by = "routeF",
            multiple = "all") 

breaks <- c(-7, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 7)
lgnd_head <- "Mean Trend\n"
trend_title <- "Mean Trend"
labls = c(paste0("< ",breaks[1]),paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),paste0("> ",breaks[length(breaks)]))
labls = paste0(labls, " %/year")
plot_map$Tplot <- cut(plot_map$trend,breaks = c(-Inf, breaks, Inf),labels = labls)


map_palette <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf",
                 "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")
names(map_palette) <- labls


scbar <- scalebar(plot_map,
                  dist = 300,
                  dist_unit = "km",
                  transform = FALSE,
                  facet.var = "parameter",
                  facet.lev = "Full with Habitat-Change",
                  location = "bottomleft",
                  st.size = 2.5,
                  box.fill = c(gray(0.5),"white"),
                  box.color = gray(0.3),
                  st.color = gray(0.5))


map <- ggplot()+
  geom_sf(data = base_strata_map,
          fill = NA,
          colour = grey(0.75))+
  geom_sf(data = plot_map,
          aes(colour = Tplot,
              size = abundance))+
  scale_size_continuous(range = c(0.05,2),
                        name = "Mean Count")+
  scale_colour_manual(values = map_palette, aesthetics = c("colour"),
                      guide = guide_legend(reverse=TRUE),
                      name = paste0(lgnd_head))+
  coord_sf(xlim = xlms,ylim = ylms)+
  guides(size = "none")+
  # scalebar(plot_map,
  #          dist = 300,
  #          dist_unit = "km",
  #          transform = FALSE,
  #          facet.var = "parameter",
  #          facet.lev = "Full with Habitat-Change",
  #          location = "bottomleft",
  #          st.size = 2.5,
  #          box.fill = c(gray(0.5),"white"),
  #          box.color = gray(0.3),
  #          st.color = gray(0.5))+
  xlab("")+
  ylab("")+
  # ggspatial::annotation_north_arrow(data = base_strata_map,
  #                                   aes(location = "tr"),
  #                                   style = north_arrow_minimal(
  #                                     line_width = 1,
  #                                     line_col = gray(0.7),
  #                                     fill = gray(0.7),
  #                                     text_col = gray(0.5),
  #                                     text_family = "",
  #                                     text_face = NULL,
  #                                     text_size = 10
  #                                   ))+
  # north(symbol = 1, location = "topright",
  #       anchor = c(x = -123,y = 61),
  #       x.min = -124,
  #       x.max = -122,
  #       y.min = 60,
  #       y.max = 62)+
  labs(title = "Trend")+
  theme_bw()+
  facet_wrap(vars(parameter))


map <- map + scbar

map


map_abund <- ggplot()+
  geom_sf(data = base_strata_map,
          fill = NA,
          colour = grey(0.75))+
  geom_sf(data = plot_map,
          aes(colour = abundance))+
  scale_colour_viridis_c(begin = 0.1, end = 0.9,
                         guide = guide_legend(reverse=TRUE),
                         name = paste0("Relative Abundance"))+
  coord_sf(xlim = xlms,ylim = ylms)+
  theme_bw()+
  xlab("")+
  ylab("")+
  # north(plot_map, symbol = 3)+
  labs(title = "Relative Abundance")+
  facet_wrap(vars(parameter))




map_se <- ggplot()+
  geom_sf(data = base_strata_map,
          fill = NA,
          colour = grey(0.75))+
  geom_sf(data = plot_map,
          aes(colour = trend_se,
              size = abundance_se))+
  scale_size_continuous(range = c(0.05,2),
                        name = "SE of Mean Count",
                        trans = "reverse")+
  scale_colour_viridis_c(aesthetics = c("colour"),
                         guide = guide_legend(reverse=TRUE),
                         name = paste0("SE of Trend"))+
  coord_sf(xlim = xlms,ylim = ylms)+
  theme_bw()+
  xlab("")+
  ylab("")+
  # north(plot_map, symbol = 3)+
  guides(size = "none")+
  labs(title = "")+
  facet_wrap(vars(parameter))
