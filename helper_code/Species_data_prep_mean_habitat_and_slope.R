## Prepare data for modeling

# This script 
# setwd("C:/Users/SmithAC/Documents/GitHub/Jefferys_etal")
library(here)
here::i_am("4 Jefferys_etal/Species_data_prep_mean_habitat_and_slope.R")
library(bbsBayes)
library(tidyverse)
library(cmdstanr)
library(patchwork)
source("functions/neighbours_define_voronoi.R") ## function to define neighbourhood relationships for spatial model comparison
source("functions/prepare-data-alt.R") ## small alteration of the bbsBayes function
## above source() call over-writes the bbsBayes prepare_data() function
source("functions/get_basemap_function.R") ## loads one of the bbsBayes strata maps
source("functions/posterior_summary_functions.R") ## functions similar to tidybayes that work on cmdstanr output



output_dir <- "output"



# load and stratify RUHU data ---------------------------------------------
strat = "bbs_usgs"
model = "slope"

strat_data <- bbsBayes::stratify(by = strat)



# for(firstYear in c(1985,2006)){
#   
#   
# lastYear <- ifelse(firstYear == 1985,2005,2021)

firstYear <- 2006
lastYear <- 2021

spp <- "_base_"

species <- "Rufous Hummingbird" 
  
  
  
species_f <- gsub(gsub(species,pattern = " ",replacement = "_",fixed = T),pattern = "'",replacement = "",fixed = T)




out_base <- paste0(species_f,spp,firstYear,"_",lastYear)
sp_data_file <- paste0("Data/",species_f,"_",firstYear,"_",lastYear,"_stan_data.RData")









# Prepare route-level model data ------------------------------------------



### this is the alternate prepare data function - modified from bbsBayes
jags_data = prepare_data(strat_data = strat_data,
                              species_to_run = species,
                              model = model,
                              #n_knots = 10,
                              min_year = firstYear,
                              max_year = lastYear,
                              min_n_routes = 1,
                         min_max_route_years = 1)



# strata map of one of the bbsBayes base maps
# helps group and set boundaries for the route-level neighbours,
## NOT directly used in the model
strata_map  <- get_basemap(strata_type = strat,
                           transform_laea = TRUE,
                           append_area_weights = FALSE)


#create list of routes and locations to ID routes that are not inside of original strata (some off-shore islands)
route_map1 = unique(data.frame(route = jags_data$route,
                              strat = jags_data$strat_name,
                              Latitude = jags_data$Latitude,
                              Longitude = jags_data$Longitude))

#create spatial object of above
route_map1 = st_as_sf(route_map1,coords = c("Longitude","Latitude"))
st_crs(route_map1) <- 4269 #NAD83 commonly used by US federal agencies

#reconcile the projections of routes and base bbs strata
route_map1 = st_transform(route_map1,crs = st_crs(strata_map))

#drops the routes geographically outside of the strata (some offshore islands) 
# and adds the strat indicator variable to link to model output
strata_map_buf <- strata_map %>% 
  filter(ST_12 %in% route_map1$strat) %>% 
  summarise() %>% 
  st_buffer(.,10000) #drops any routes with start-points > 10 km outside of strata boundaries
realized_routes <- route_map1 %>% 
  st_join(.,strata_map_buf,
          join = st_within,
          left = FALSE) 



# reorganizes data after routes were dropped outside of strata
new_data <- data.frame(strat_name = jags_data$strat_name,
                       strat = jags_data$strat,
                       route = jags_data$route,
                       strat = jags_data$strat_name,
                       Latitude = jags_data$Latitude,
                       Longitude = jags_data$Longitude,
                       count = jags_data$count,
                       year = jags_data$year,
                       firstyr = jags_data$firstyr,
                       ObsN = jags_data$ObsN,
                       r_year = jags_data$r_year) %>% 
  filter(route %in% realized_routes$route)





# add route-names
rt_names_counts <- strat_data$route_strat %>% 
  select(statenum,Route,RouteName,rt.uni,State) %>% 
  distinct() %>% 
  filter(rt.uni %in% unique(new_data$route)) 

# Link covariates and counts ----------------------------------------------

# Load hsi data -------------------------------------------------------

hsi <- read.csv("data (original)/MeanHSI_BBSRoute_yearly.csv")



# clean route names to separate names and numbers in canadian routes
hsi <- hsi %>% 
  mutate(RTENO_2 = ifelse(str_detect(RTENAME,"[[:digit:]]"),
                          str_extract_all(RTENAME,"[[:digit:]]{2}-[[:digit:]]{3}"),
                          RTENO),
         RTENO_2 = str_replace(RTENO_2,"-",
                               ""),
         RouteName = ifelse(str_detect(RTENAME,"[[:digit:]]"),
                            str_replace(RTENAME,"[[:digit:]]{2}-[[:digit:]]{3}[[:space:]]",""),
                            RTENAME),
         rt.uni = paste0(as.integer(str_sub(RTENO_2,1,2)),
                         "-",
                         as.integer(str_sub(RTENO_2,3,5))),
         year = as.integer(str_sub(date,1,4))) 


# hab_vis <- ggplot(data = hsi,
#                   aes(x = year,y = meanHSI,
#                       group = rt.uni))+
#   geom_line(alpha = 0.2)
# 
# hab_vis

#

### joint the list of routes with BBS survey data and the habitat info
join_test1 <- rt_names_counts %>% 
  filter(statenum != 3) %>% 
  left_join(.,hsi,
            by = c("rt.uni","RouteName"),
            multiple = "all")

### Alaska route numbers are not in the file, so need to join
### separately.
AK <- rt_names_counts %>% 
  filter(statenum == 3) %>% 
  left_join(.,hsi,
            by = "RouteName",
            multiple = "all") %>% 
  select(-rt.uni.y) %>% 
  rename(rt.uni = rt.uni.x) 
### combine AK and rest of route names and numbers
join_test <- bind_rows(AK,join_test1)


# there are some routes for which route-GIS data were not available
missing_habitat <- join_test %>% 
  filter(is.na(meanHSI))

write.csv(missing_habitat,
          paste0("routes_with_RUHU_data_missing_habitat",firstYear,".csv"))


### drop the routes with no HSI information
join_table <- join_test %>% 
  filter(!is.na(meanHSI))

hab_full <- join_table
### hab_full is a full list of annual hsi measures

# hab_vis <- ggplot(data = hab_full,
#                   aes(x = year,y = meanHSI,
#                       group = RouteName))+
#   geom_line(alpha = 0.2)+
#   geom_point(alpha = 0.1)+
#   geom_smooth(method = "lm",se = FALSE,
#               alpha = 0.2,
#               linewidth = 0.2)+
#   facet_wrap(vars(State))
# 
# hab_vis



hab_full <- hab_full %>% 
  rename(route = rt.uni) %>% 
  select(route,year,meanHSI) %>% 
  filter(year >= firstYear)





# Drop count data with no habitat -----------------------------------------

new_data <- new_data %>% 
  inner_join(.,hab_full,
             by = c("route",
                    "r_year" = "year"))

# hab_vis <- ggplot(data = new_data,
#                   aes(x = r_year,y = meanHSI,
#                       group = route))+
#   geom_line(alpha = 0.2)+
#   geom_point(alpha = 0.1)+
#   geom_smooth(method = "lm",se = FALSE,
#               alpha = 0.2,
#               linewidth = 0.2)+
#   facet_wrap(vars(strat_name))
# 
# hab_vis




# Spatial neighbours set up --------------------
strata_list <- data.frame(ST_12 = unique(new_data$strat_name),
                          strat = unique(new_data$strat))


realized_strata_map <- strata_map %>%
  filter(ST_12 %in% strata_list$ST_12)

new_data$routeF <- as.integer(factor((new_data$route))) #main route-level integer index

#create a data frame of each unique route in the species-specific dataset
route_map = unique(data.frame(route = new_data$route,
                              routeF = new_data$routeF,
                              strat = new_data$strat_name,
                              Latitude = new_data$Latitude,
                              Longitude = new_data$Longitude))


# reconcile duplicate spatial locations -----------------------------------
# adhoc way of separating different routes with the same starting coordinates
# this shifts the starting coordinates of the duplicates by ~1.5km to the North East 
# ensures that the duplicates have a unique spatial location, but remain very close to
# their original location and retain reasonable neighbourhood relationships
# these duplicates happen when a "new" route is established (i.e., a route is re-named) 
# because some large proportion
# of the route is changed, but the start-point remains the same
dups = which(duplicated(route_map[,c("Latitude","Longitude")]))
while(length(dups) > 0){
  route_map[dups,"Latitude"] <- route_map[dups,"Latitude"]+0.01 #=0.01 decimal degrees ~ 1km
  route_map[dups,"Longitude"] <- route_map[dups,"Longitude"]+0.01 #=0.01 decimal degrees ~ 1km
  dups = which(duplicated(route_map[,c("Latitude","Longitude")]))
  
}
dups = which(duplicated(route_map[,c("Latitude","Longitude")])) 
if(length(dups) > 0){stop(paste(spec,"At least one duplicate route remains"))}



#create spatial object from route_map dataframe
route_map = st_as_sf(route_map,coords = c("Longitude","Latitude"))
st_crs(route_map) <- 4269 #NAD83 commonly used by US federal agencies


#reproject the routes spatial object ot match the strata-map in equal area projection
route_map = st_transform(route_map,crs = st_crs(strata_map))

## custom function that returns the adjacency data necessary for the stan model
## also exports maps and saved data objects to plot_dir

## this function generates a Voronoi tesselation with each route start-point
## at teh centre of the voronoi polygons.
## it then clips the tesselated surface using an intersection of the strata map 
## and a concave polygon surrounding the route start locations
## this clipping ensures that routes at the edges of the species range are not able
## to be considered neighbours with other very-distant edge-routes and retains
## the complex shape of many species distributions (e.g., species with boreal and 
## mountainous distributions, such as Pileated Woodpecker)
## it also forces a completely connected network of routes. So if some portions of the 
## species range are separated by gaps (e.g., intervening BBS strata without any routes)
## it forces a neighbourhood connection between the closest pair of routes that would
## bridge the gap (and any other routes within 10% of the distance of the closest pair)
car_stan_dat <- neighbours_define_voronoi(real_point_map = route_map,
                                  species = species,
                                  strat_indicator = "routeF",
                                  strata_map = realized_strata_map,
                                  concavity = 1)#concavity argument from concaveman()
## a relative measure of concavity. 1 results in a relatively detailed shape, Infinity results in a convex hull.

## see the map of connections
print(car_stan_dat$map)


pdf(paste0("data/maps/route_map_",firstYear,"-",lastYear,"_",species_f,".pdf"))
print(car_stan_dat$map)
dev.off()


# habitat predictor data -------------------------------------------

route_df <- route_map %>% 
  sf::st_drop_geometry(route_map)


sl <- function(y,x){
  t <- lm(y~x)
  cc <- coefficients(t)[["x"]]
  return(cc)
}
slope_full <- NULL
 
  
  hab_route_mean_t <- hab_full %>% 
    #filter(year < (firstYear + 10)) %>%  #first 10 years
    group_by(route) %>% 
    summarise(mean_hsi_f6 = mean(meanHSI,na.rm = T))

  hab_route_slope_t <- hab_full %>% 
    arrange(year,route) %>% 
    group_by(route) %>% 
    summarise(slope_hsi = sl(meanHSI,year))

  hab_tmp <- inner_join(hab_route_mean_t,hab_route_slope_t,
                        by = "route") %>% 
    mutate(startYear = firstYear)
  
  slope_full <- bind_rows(slope_full,hab_tmp)
    


slope_full <- slope_full %>% 
  left_join(.,route_df,
            by = "route") %>% 
  arrange(routeF)
  

### Build the data list required for Stan


stan_data <- list()
stan_data[["count"]] <- new_data$count
stan_data[["ncounts"]] <- length(new_data$count)
stan_data[["strat"]] <- new_data$strat
stan_data[["route"]] <- new_data$routeF
stan_data[["year"]] <- new_data$year
stan_data[["firstyr"]] <- new_data$firstyr
stan_data[["fixedyear"]] <- jags_data$fixedyear # mis-year of time series


stan_data[["nyears"]] <- max(new_data$year)
stan_data[["observer"]] <- as.integer(factor((new_data$ObsN)))
stan_data[["nobservers"]] <- max(stan_data$observer)



stan_data[["N_edges"]] <- car_stan_dat$N_edges
stan_data[["node1"]] <- car_stan_dat$node1
stan_data[["node2"]] <- car_stan_dat$node2
stan_data[["nroutes"]] <- max(stan_data$route)
stan_data[["route_habitat"]] <- as.numeric(slope_full$mean_hsi_f6 - mean(slope_full$mean_hsi_f6))
stan_data[["route_habitat_slope"]] <- 100*(slope_full$slope_hsi - mean(slope_full$slope_hsi))
 

if(car_stan_dat$N != stan_data[["nroutes"]]){stop("Some routes are missing from adjacency matrix")}

dist_matrix_km <- dist_matrix(route_map,
                              strat_indicator = "routeF")
save(list = c("stan_data",
              "new_data",
              "route_map",
              "realized_strata_map",
              "car_stan_dat",
              "dist_matrix_km",
              "hab_full",
              "slope_full"),
     file = sp_data_file)


mean_obs <- new_data %>% 
  filter(r_year < firstYear+15) %>% 
  group_by(route) %>% 
  summarise(mean_obs = mean(count,na.rm = TRUE),
            .groups = "keep")


  
rt_sums <- slope_full %>% 
  left_join(.,mean_obs,
            by = "route") 

map_df <- route_map %>% 
  left_join(rt_sums,by = c("route"))


t1 <- ggplot(data = map_df)+
  geom_sf(aes(colour = slope_hsi,
              size = mean_obs))+
  colorspace::scale_colour_continuous_diverging(mid = 0,
                                                rev = TRUE)+
  labs(title = paste0("HSI trend ",firstYear,"-",lastYear,"\n by mean counts in ",firstYear))+
  theme_bw()



pdf(paste0("Figures/HSI_trend_by_obs_counts_",firstYear,".pdf"),
    width = 8,height = 8)
print(t1)
dev.off()


# }


