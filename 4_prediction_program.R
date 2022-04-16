################################################################################
# ABOUT THIS SCRIPT
################################################################################
# AUTHOR, DATE: 
# Magali Blanco, 1/26/2022

# OBJECTIVE
# This program takes in a dataset with locations and their respective geographic covariates, and it returns a dataset with annual average air pollution predictions for UFP, BC, NO2, PM2.5 and CO2 at those locations.
# Predictions are based on the data collected from a year-long air quality monitoring campaign in the greater Seattle area.
# The campaign had 309 stop locations, each with approximately 29 repeat visits.

# INPUTS
# Three inputs are required: 
#   1. the file path for the geocovariate dataset where predictions are desired
#   2. the repository where prediction files should be saved
#   3. the desired prediction file format, either csv or rda 

# OUTPUTS
# The output is a dataset with annual average air pollution predictions (UFP, BC, NO2, PM2.5 and CO2) for the locations with geographic covariates. 
# Plots and maps to summarize the resulting predictions are also included.

# ERROR MESSAGES
# Error messages occur if:
#   1. three arguments are not included when running this program
#   2. there are locations with missing covariate values or missing covariates altogether that are required for the prediction models to run

# EXAMPLE OF HOW TO USE THIS PROGRAM
# in a terminal open to the R program project directory, type: rscript coding_example.R <covariate_file_path> <prediction_directory> <prediction_file_format>
# rscript 4_prediction_program.R Data/Original/Geocovariates/dr0311_grid_covars.rda output/Predictions/grid csv 

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
pacman::p_load(tidyverse, ggpubr, pls, gstat, sf, ggspatial, tictoc, tools)

# report how long script takes to run
tic()

# ensure reproducibility 
set.seed(1)

###########################################################################################
# TAKE IN USER ARGUMENTS 
###########################################################################################
#allow R to take input from the command line
user_arguments <- commandArgs(trailingOnly = TRUE)
# #test 
# user_arguments <-c("Data/Original/Geocovariates/dr0311_grid_covars.rda", "Data/Output/Predictions/grid_test", "csv")

if (length(user_arguments) !=3) {
  print("Usage error. Enter: 1. the location of the covariate dataset for which you would like predictions, 2. where the prediction outputs should be saved, and 3. the desired prediction file fomat (csv or rda). Usage:")
  print("rscript coding_example.R <covariate_file_path> <prediction_directory> <prediction_file_format>")
  stop()
}

# new covariate file
covariate_file_path <- user_arguments[1]
cov_ext <- file_ext(covariate_file_path)

#where predictions should be saved
prediction_directory <- user_arguments[2]
## create the directory if it does not already exists
if(!dir.exists(prediction_directory)) {dir.create(prediction_directory)}

# the prediction ffile format (e.g., 'rda')
prediction_file_format <- tolower(user_arguments[3])

###########################################################################################
# UPLOAD THE NEW DATASET WHERE PREDICTIONS ARE DESIRED
###########################################################################################
if(cov_ext == "rda") {dt0 <- readRDS(covariate_file_path)}
if(cov_ext == "csv") {dt0 <- read_csv(covariate_file_path)}

if(!cov_ext %in% c("csv", "rda")) {stop("Error. Covariate file must be a CSV or RDA file")}

###########################################################################################
# UPLOAD MODELING DATA
###########################################################################################
# desired prediction averaging approach. The mean of winsorized medians approach takes multiple median concentration readings at each monitoring location (n=309 sites x ~29 visits/site), trims extreme observations, and calculates an overall average.
keep_averages <- c("mean_of_win_medians")

# full dataset that will be used to build a prediction model. It contains 309 sites with annual average estimates (Y) and geographic covariates (X)
annual <- readRDS(file.path("Data", "Output", "annual.rda")) %>%
  #only want this averaging approach
  filter(annual %in% keep_averages)

# the covariate names that will be used in the model
cov_names <- select(annual, log_m_to_a1:last_col()) %>%
  names()

# load a spatial file of the original monitoring area to assess spatial extrapolation later
monitoring_area <- readRDS(file.path("Data", "Output", "GIS", "monitoring_land_zero_water_shp.rda" #"monitoring_land_shp.rda"
                                     ))  
lat_long_crs <- 4326

###########################################################################################
# Universal Kriging - Partial Least Squares Model function

# desired PLS components to use (from a different script): 3
load(file.path("Data", "Output", "Objects", "pls_comp_n.rda"))

#prediction model
uk_pls <- readRDS(file.path("Data", "Output", "Predictions", "uk_pls_model.rda"))


###########################################################################################
# GENERATE NEW COVARIATES FOR THE DATASET
###########################################################################################
# created some new proximity variables  
# log transform land proximity variables (e.g., distance to roadways)

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
  # for the NO2 covariate, use the average levels from several available years
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
  dt %>%
    st_as_sf(coords = c('longitude', 'latitude'), crs= lat_long_crs) %>%
    st_intersects(., monitoring_area, sparse = F) %>%
    apply(., 1, any)
)

###########################################################################################
# QC CHECKS
###########################################################################################
# 1. ensure that the new dataset has all the necessary geograpahic covariates for modeling, otherwise give an error message with the missing covariates
has_all_covariates <- all(cov_names %in% names(dt))

if(has_all_covariates==FALSE) {
  missing_cov <- cov_names[!cov_names %in% names(dt)]
  error_msg <- paste("The following covariates are needed but missing from the dataset:", paste(missing_cov, collapse = ", "), ". Please fix this before continuing.")
  stop(error_msg)
}

# 2. check that there are no missing covariate values for any location
has_missing_values <- sapply(dt[cov_names], function(x) any(is.na(x) )) %>% 
  as.data.frame() %>% rownames_to_column() %>%
  rename(true_or_false = '.')  

if(any(has_missing_values$true_or_false) == TRUE) {
  covariates_with_missingness <- filter(has_missing_values, true_or_false==TRUE) %>%
    pull(rowname)
  
  error_msg <- paste("The following covariates have missing values:", paste(covariates_with_missingness, collapse = ", "), ". Please fix this before continuing.")
  
  stop(error_msg)
}

# 3. print a 'pass' message if all of the covariates are present and there are no locations with missing values
if(has_all_covariates ==TRUE & any(has_missing_values$.) == FALSE) {print("Covariate checks passed.")} 

###########################################################################################
# PREDICT AT NEW DATASET
###########################################################################################
print("Generating predictions...")

# for each pollutant ('variable') and annual average estimate ('annual'), build prediction models using the monitorign data and predict at the new dataset locations
new_predictions0 <- lapply(group_split(annual, variable, annual), function(x) {
  temp <- dt %>%
    mutate(
      variable = first(x$variable),
      annual = first(x$annual)
    ) %>%
    uk_pls(new_data = ., modeling_data = x)
}) %>%
  bind_rows()  

# save the location and prediction information
new_predictions <- new_predictions0 %>%
  select(contains(c("_id", "_key", "msa")), longitude, latitude, in_monitoring_area, variable, annual, prediction)

###########################################################################################
# SAVE THE PREDICTIONS
###########################################################################################
prediction_file_name <- file.path(prediction_directory, paste0("predictions.", prediction_file_format))

if(prediction_file_format == "csv") {write.csv(new_predictions, prediction_file_name, row.names = F)}
if(prediction_file_format == "rda") {saveRDS(new_predictions, prediction_file_name)}

print(paste0("Predictions saved: ", prediction_file_name))

###########################################################################################
# GENERATE SUMMARY FIGURES AND MAPS OF THE NEW PREDICTIONS
###########################################################################################
print("Generating prediction summary figures and maps...")

# 1. prediction histograms 
ggplot(data=new_predictions, aes(x=prediction)) + 
  facet_wrap(~variable, scales = "free") + 
  geom_histogram(bins=30) +
  labs(title = "Prediction Histograms")

ggsave(file.path(prediction_directory, "prediction_histograms.png"), height = 10, width = 8)

# 2. prediction maps at all of the locations
p <- list()

for (i in unique(new_predictions$variable)) {
  
  df <- filter(new_predictions, variable == i)
  
  p[[i]] <- ggplot() +
    geom_sf(data=monitoring_area)  +
    geom_point(data= df, aes(x=longitude, y=latitude, col=prediction), alpha=0.6) + 
    facet_wrap(~variable) +   
    scale_colour_gradient(name = "Conc", low = "yellow", high = "red") +
    # add scales & N arrow 
    annotation_scale(location = "tr", unit_category ="imperial", pad_y = unit(0.55, "cm")) +
    annotation_north_arrow(location = "tr", pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering) +
    theme_bw() +
    coord_sf(expand = F) 
}

ggarrange(plotlist = p) %>% 
  annotate_figure(top = "UK-PLS predictions at all locations")

ggsave(file.path(prediction_directory, "all_predictions.png"), height = 16, width = 16)

# 3. prediction maps only at locations in monitoring area
p <- list()

for (i in unique(new_predictions$variable)) {
  
  df <- filter(new_predictions,
               in_monitoring_area == TRUE,
               variable == i)
  
  p[[i]] <- ggplot() +
    geom_sf(data=monitoring_area)  +
    geom_point(data= df, aes(x=longitude, y=latitude, col=prediction), alpha=0.6) + 
    facet_wrap(~variable) + 
    scale_colour_gradient(name = "Conc", low = "yellow", high = "red") +
    # add scales & N arrow 
    annotation_scale(location = "tr", unit_category ="imperial", pad_y = unit(0.55, "cm")) +
    annotation_north_arrow(location = "tr", pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering) +
    theme_bw() +
    coord_sf(expand = F)
}

ggarrange(plotlist = p) %>% 
  annotate_figure(top = "UK-PLS predictions at locations in monitoring area")

ggsave(file.path(prediction_directory, "monitoring_predictions.png"), height = 16, width = 16)

print(paste0("Figures and maps saved to the following directory: ", prediction_directory))

###########################################################################################
#print the script run duration
toc()

print("Program done.")

###########################################################################################
