# Magali Blanco
# 4/17/22

# Script purpose:
## 1. this script cleans the cohort covariates to only include locations within the monitoring or study area
## 2. it includes code to run in the terminal to make cohort predictions

################################################################################
# SETUP
################################################################################
# Clear workspace of all objects and unload all extra (non-base) packages
rm(list = ls(all = TRUE))
if (!is.null(sessionInfo()$otherPkgs)) {
  res <- suppressWarnings(
    lapply(paste('package:', names(sessionInfo()$otherPkgs), sep=""),
           detach, character.only=TRUE, unload=TRUE, force=TRUE))
}

# load the required libraries for: plotting, modeling, spatial features, script timing
if (!require("pacman")) {install.packages("pacman")}
pacman::p_load(tidyverse, sf, ggspatial, readr)

# ensure reproducibility 
set.seed(1)

################################################################################
# Data
################################################################################
dt <- read_csv("Data/Original/Geocovariates/dr0357_cohort_covar_20220404.csv")

# load a spatial file of the original monitoring area to assess spatial extrapolation later
monitoring_area <- readRDS(file.path("Data", "Output", "GIS", "monitoring_land_zero_water_shp.rda" #"monitoring_land_shp.rda"
))  
lat_long_crs <- 4326


###########################################################################################
# ADD LOCATION INDICATORS
###########################################################################################
# indicators of whether or not the prediction locations are in the monitoring area
in_monitoring_area <- suppressMessages(
  dt %>%
    st_as_sf(coords = c('longitude', 'latitude'), crs= lat_long_crs) %>%
    st_intersects(., monitoring_area, sparse = F) %>%
    apply(., 1, any)
)

###########################################################################################
# ONLY KEEP LOCATIONS IN THE MONITORING REGION
###########################################################################################
dt2 <- filter(dt, in_monitoring_area)


# # check that there are no missing pop or bus values anymore
# dt2 %>% 
#   mutate(missing = is.na(bus_s00050)) %>%
#   ggplot(aes(x=longitude, y=latitude,col=missing)) + 
#   geom_point(alpha=0.05)
#   

###########################################################################################
# SAVE DATA
###########################################################################################
write.csv(dt2, file.path("Data", "Output", "dr0357_cohort_covar_20220404_in_monitoring_area.csv"), row.names = F)



###########################################################################################
# TEST
###########################################################################################
# # check missing data
# test <- read_csv(file.path("Data", "Output", "Predictions", "Cohort", "predictions.csv"))
# 
# # missing prediction locations are always outside the monitoring region
# table(monitoring_area = test$in_monitoring_area, missing_prediction = is.na(test$prediction))

