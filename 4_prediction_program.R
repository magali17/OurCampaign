
#EXAMPLE OF HOW TO USE THIS PROGRAM
# in a terminal open to the R project directory type:
# rscript 4_prediction_program.R Data/Original/Geocovariates/dr0311_grid_covars.rda Data/Output/Predictions/grid rda 

################################################################################
# OBJECTIVE
# This program takes in a dataset with locations and their respective geographic covariates, and it returns a dataset with annual average air pollution predictions for PNC, BC, NO2, PM2.5 and CO2 at those locations.

# INPUTS
# Three inputs are required: 
#   1. the file path for the geocovariate dataset where predictions are desired
#   2. the repository where prediction files should be saved
#   3. the desired prediction file format, either csv or rda 

# OUTPUTS
# The output is a dataset with annual average air pollution predictions (PNC, BC, NO2, PM2.5 and CO2) for the locations with geograhic covariates. 
# Plots and maps to summarize the resulting predictions are also included.
# The default annual averages modeled are the mean of windzorized medians

# ERROR MESSAGES
# Error messages occur if:
#   1. three arguments are not included
#   2. there are locations with missing covariates or missing covariates alltogether that are required for the prediction models to run

################################################################################

# Clear workspace of all objects and unload all extra (non-base) packages
rm(list = ls(all = TRUE))
if (!is.null(sessionInfo()$otherPkgs)) {
  res <- suppressWarnings(
    lapply(paste('package:', names(sessionInfo()$otherPkgs), sep=""),
           detach, character.only=TRUE, unload=TRUE, force=TRUE))
}

pacman::p_load(tidyverse, ggpubr,
               pls, gstat,  
               sf, #parallel,
               tictoc
)

#report how long script takes to run
tic()

set.seed(1)

#functions
source("0_Functions.R")

#use_cores <- 4

###########################################################################################
# USER ARGUMENTS 
###########################################################################################
#allow R to take input from the command line
user_arguments <- commandArgs(trailingOnly = TRUE)
#test 
#user_arguments <-c("Data/Original/Geocovariates/dr0311_grid_covars.rda", "Data/Output/Predictions/grid_test", "csv")

if (length(user_arguments) !=3) {
  print("Usage error. Enter: 1. the location of the covariate dataset for which you would like predictions, 2. where the prediction outputs should be saved, and 3. the desired prediction file fomat (csv or rda). Usage:")
  print("rscript prediction_program.R <covariate_file_path> <prediction_directory> <prediction_file_format>")
  stop()
}


###########################################################################################
# FILE DETAILS
###########################################################################################
# use file.choose() to get full file path 

# new covariate file
covariate_file_path <- user_arguments[1]

#where predictions should be saved
prediction_directory <- user_arguments[2]
## create the directory if it does not already exists
if(!dir.exists(prediction_directory)) {dir.create(prediction_directory)}

# the prediction fiel format (e.g., 'rda')
prediction_file_format <- tolower(user_arguments[3])

# desired predictions averages
keep_averages <- c("mean_of_win_medians")
#save(keep_averages, file = file.path("Data", "Output", "Objects", "keep_averages.rda"))

###########################################################################################
# UPLOAD NEW DATASET
###########################################################################################
dt0 <- readRDS(covariate_file_path)

###########################################################################################
# UPLOAD MODELING DATA
###########################################################################################
# mapping variables
lat_long_crs <- 4326

# full dataset: 309 site estimates and geocovariates
annual <- readRDS(file.path("Data", "Output", "annual.rda")) %>%
  #only want these averages
  filter(annual %in% keep_averages)

#prediction model
uk_pls <- readRDS(file.path("Data", "Output", "Predictions", "uk_pls_model.rda"))
#pls_comp_n <-2
# use 3 PLS components, or  could be different 
load(file.path("Data", "Output", "Objects", "pls_comp_n.rda"))

# covariate names
cov_names <- annual %>%
  select(log_m_to_a1:last_col()) %>%
  names()

#shapefiles
monitoring_area <- readRDS(file.path("Data", "Output", "GIS", "monitoring_land_shp.rda"))  

###########################################################################################
# NEW DATASET COVARIATES
###########################################################################################
# FN TO ADD NEW COVARIATES TO NEW DATASET

# * transform and generate new variables for the monitoring and all other relevant data (agency sites, cohort sites), .  
#    + log transformed proximities variables  
#       -	TRAP tends to exponentially decays with increasing distance from the source    
#    + created some new proximity variables (e.g. distance to a123)

combine_a23_ll <- function(df) {
  #find buffers for a2-3 length variables
  buffers <- str_subset(names(df), "ll_a[2:3]") %>% str_extract("s[0:9].*")
  
  #for each buffer, calculate sum of a2+a3 length
  for (i in seq_along(buffers)) {
    old_vars <- paste0(c("ll_a2_", "ll_a3_"), buffers[i])
    new_var <- paste0("ll_a23_", buffers[i])
    
    df[new_var] <- apply(df[old_vars], 1, sum)
  }
  return(df)
}

generate_new_vars <- function(df) {
  no2_behr_vars <- c("no2_behr_2005","no2_behr_2006", "no2_behr_2007")
  
  df <- df %>%
    rowwise() %>%
    mutate(m_to_a123 = min(m_to_a1, m_to_a2, m_to_a3),
           m_to_a23 = min(m_to_a2, m_to_a3),
           
           no2_behr = mean(!!as.symbol(no2_behr_vars))
    ) %>%
    ungroup() %>%
    #make min distance 1 m before log transforming
    mutate_at(vars(starts_with("m_to_")), ~ifelse(.==0, 1, .) %>% log(.)) %>%
    rename_at(vars(starts_with("m_to_")), ~gsub("m_to_", "log_m_to_", .)) %>%
    # calculate sum of a2 and a3 roads in each buffer
    combine_a23_ll()
}

###########################################################################################

dt <- generate_new_vars(dt0)

###########################################################################################
# ADD LOCATION INDICATORS
###########################################################################################

# add indicators of whether or not the prediction locations are in the monitoring area
dt$in_monitoring_area <- suppressMessages(
  dt  %>%
    st_as_sf(coords = c('longitude', 'latitude'), crs=lat_long_crs) %>%
    st_intersects(., monitoring_area, sparse = F) %>%
    apply(., 1, any)
  )

###########################################################################################
# QC CHECKS
###########################################################################################
# check that new dataset has all the geocovariates necessary for modeling

has_all_covariates <- all(cov_names %in% names(dt))

if(has_all_covariates==FALSE) {
  missing_cov <- cov_names[!cov_names %in% names(dt)]
  error_msg <- paste("The following covariates are needed but missing from the dataset:", paste(missing_cov, collapse = ", "), ". Please fix this before continuing.")
  stop(error_msg)
}

#check that there are no missing values
has_missing_values <- sapply(dt[cov_names], function(x) any(is.na(x) )) %>% 
  as.data.frame() %>% rownames_to_column() 

if(any(has_missing_values$.) ==TRUE) {
  covariates_with_missingness <- has_missing_values %>%
    filter(.==TRUE) %>%
    pull(rowname)
  
  error_msg <- paste("The following covariates have missing values:", paste(covariates_with_missingness, collapse = ", "), ". Please fix this before continuing.")
  
  stop(error_msg)
}
# print a message if all of the covariates are present and there are no locations with missing values
if(has_all_covariates ==TRUE & any(has_missing_values$.) == FALSE) {print("Covariate checks passed.")} 

###########################################################################################
# PREDICT AT NEW DATASET
###########################################################################################
print("Generating predictions...")

# build prediction models for each pollutant ('variable') and annual average estimate; predict at new locations
new_predictions0 <- lapply(group_split(annual, variable, annual), function(x) {
    temp <- dt %>%
      mutate(
        variable = first(x$variable),
        annual = first(x$annual)
        ) %>%
      uk_pls(new_data = ., modeling_data = x)
    }#, mc.cores = use_cores
    ) %>%
  bind_rows()  
  

new_predictions <- new_predictions0 %>%
  select(contains(c("_id", "_key", "msa")), longitude, latitude, in_monitoring_area, variable, annual, prediction)

###########################################################################################
# CHECK THAT ALL LOCATIONS HAVE PREDICTIONS
###########################################################################################

# print("missing predictions")
# 
# new_predictions %>%
#   group_by(variable, annual) %>%
#   summarize(missing_predictions = sum(is.na(prediction)))
# 
# print("unique locations with predictions") 
# new_predictions %>%
#   distinct(longitude, latitude) %>%
#   nrow()

###########################################################################################
# FIGURES/MAPS OF PREDICTIONS
###########################################################################################
# histogram of predictions
ggplot(data=new_predictions, aes(x=prediction)) + 
  facet_wrap(~variable, scales = "free") + 
  geom_histogram(bins=30) +
  labs(title = "Prediction Histograms")

ggsave(file.path(prediction_directory, "prediction_histograms.png"),width = 8, height = 10 )


# prediction maps at all of the locations
p <- list()

for (i in unique(new_predictions$variable)) {

  df <- filter(new_predictions, variable == i)
  
  p[[i]] <- ggplot() +
    geom_sf(data=monitoring_area)  +
    geom_point(data= df, aes(x=longitude, y=latitude, col=prediction), alpha=0.3) + 
    facet_wrap(~variable) +   
    theme_void()
}

ggarrange(plotlist = p) %>% 
  annotate_figure(top = "UK-PLS predictions for all locations")

ggsave(file.path(prediction_directory, "all_predictions.png"), height = 11, width = 8)

# prediction maps only at locations in monitoring area
p <- list()

for (i in unique(new_predictions$variable)) {

  df <- filter(new_predictions,
               in_monitoring_area == TRUE,
               variable == i)
  
  p[[i]] <- ggplot() +
    geom_sf(data=monitoring_area)  +
    geom_point(data= df, aes(x=longitude, y=latitude, col=prediction), alpha=0.3) + 
    facet_wrap(~variable) + 
    theme_void()
}

ggarrange(plotlist = p) %>% 
  annotate_figure(top = "UK-PLS predictions for locations in monitoring area")

ggsave(file.path(prediction_directory, "monitoring_predictions.png"), height = 11, width = 8)

###########################################################################################
# SAVE PREDICTIONS
###########################################################################################
prediction_file_name <- file.path(prediction_directory, paste0("predictions.", prediction_file_format))
                             
if(prediction_file_format == "csv") {write.csv(new_predictions, prediction_file_name, row.names = F)}
if(prediction_file_format == "rda") {saveRDS(new_predictions, prediction_file_name)}
 
###########################################################################################

print("Done.")

#print script run duration
toc()

###########################################################################################
 