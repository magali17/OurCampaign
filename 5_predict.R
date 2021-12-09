# Clear workspace of all objects and unload all extra (non-base) packages
rm(list = ls(all = TRUE))
if (!is.null(sessionInfo()$otherPkgs)) {
  res <- suppressWarnings(
    lapply(paste('package:', names(sessionInfo()$otherPkgs), sep=""),
           detach, character.only=TRUE, unload=TRUE, force=TRUE))
}

pacman::p_load(#knitr, kableExtra,
               tidyverse, ggpubr,
               pls, gstat, #variogram()
               sf,
               parallel
               )

set.seed(1)

#functions
source("0_Functions.R")

use_cores <- 4

###########################################################################################
# UPLOAD NEW DATASET
###########################################################################################

#new_file_path <- file.path("..", "..", "Kaya", "block10_intpts_wa.rda") 
new_file_path <- file.path("..", "..", "Kaya", "old", "dr0342_block_covars.rda") 

if(file.exists(new_file_path)) {
  dt0 <- readRDS(new_file_path)
  } else {
    dt0 <- read.csv(gsub(".rda", ".txt|.csv", new_file_path))
    saveRDS(dt0, file.path(new_file_path))
    }

#desired predictions
keep_averages <- c("median_of_means")
#keep_pollutants <- str_subset(unique(annual$variable), "" )

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
pls_comp_n <-2

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
  
  df <- df %>%
    rowwise() %>%
    mutate(m_to_a123 = min(m_to_a1, m_to_a2, m_to_a3),
           m_to_a23 = min(m_to_a2, m_to_a3)) %>%
    ungroup() %>%
    #make min distance 1 m before log transforming
    mutate_at(vars(starts_with("m_to_")), ~ifelse(.==0, 1, .) %>% log(.)) %>%
    rename_at(vars(starts_with("m_to_")), ~gsub("m_to_", "log_m_to_", .)) %>%
    # calculate sum of a2 and a3 roads in each buffer
    combine_a23_ll()
  }

###########################################################################################

dt <- dt0 %>%
  # NAs in bus buffers are supposed to be 0s
  mutate_at(vars(starts_with("bus_s")), ~ifelse(is.na(.), 0, .)) %>%
  generate_new_vars()


###########################################################################################
# ADD LOCATION INDICATORS
###########################################################################################

#add indicators of where prediction locations are
dt$in_monitoring_area <- dt  %>% 
  st_as_sf(coords = c('longitude', 'latitude'), crs=lat_long_crs) %>%
  st_intersects(., monitoring_area, sparse = F ) %>% 
  apply(., 1, any) #%>% st_drop_geometry()



###########################################################################################
# QC CHECKS
###########################################################################################
# check that new dataset has all the geocovariates necessary for modeling

has_all_covariates <- all(cov_names %in% names(dt))

if(has_all_covariates==FALSE) {
  missing_cov <- cov_names[!cov_names %in% names(dt)]
  error_msg <- paste("the following covariates are needed but missing from the dataset:", paste(missing_cov, collapse = ", ") )
  stop(error_msg)
}


#check that there are no missing values

has_missing_values <- sapply(dt[cov_names], function(x) any(is.na(x) )) %>% 
  as.data.frame() %>% rownames_to_column() 

if(any(has_missing_values$.) ==TRUE) {
  covariates_with_missingness <- has_missing_values %>%
    filter(.==TRUE) %>%
    pull(rowname)
  
  error_msg <- paste("the following covariates have missing values:", paste(covariates_with_missingness, collapse = ", ") )
  
  stop(error_msg)
  }

## --> add, if passes both, print("covariate checks passed"); else, print("___ not passed")


# ggplot() + geom_sf(data=monitoring_area)  +
#   geom_point(data= filter(dt, in_monitoring_area ==TRUE), # %>% mutate(missing_bus_variables = is.na(bus_s00750)), 
#              aes(x=longitude, y=latitude, col=bus_s00050), alpha=0.3) 

###########################################################################################
# PREDICT AT NEW DATASET
###########################################################################################

# x = group_split(annual, variable, annual)[[1]]
# build prediction models for each pollutant and annual average estimate; predict at new locations
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

### --> why so many missing values???

new_predictions %>%
group_by(variable, annual) %>%
  summarize(missing_predictions = sum(is.na(prediction)))

 
###########################################################################################
# FIGURES/MAPS OF PREDICTIONS
###########################################################################################
#histogram of predictions
ggplot(data=new_predictions, aes(x=prediction)) + 
  facet_wrap(~variable, scales = "free") + 
  geom_histogram() +
  labs(title = "Prediction Histograms")

ggsave(file.path("..", "..", "Kaya", "Block Predictions", "prediction_histograms.png"), )


# all locations
ggplot() +
  geom_point(data= filter(new_predictions, 
                          grepl("noscreen", variable), #!is.na(prediction)
                          ), 
             aes(x=longitude, y=latitude, col=prediction), alpha=0.3) + 
  geom_sf(data=monitoring_area, aes(fill="monitoring"), alpha=0)  +
  labs(title = "example: pnc_noscreen UK-PLS predictions for all locations",
       fill = "")

ggsave(file.path("..", "..", "Kaya", "Block Predictions", "all_predictions.png"), )

# only locations in monitoring area
p <- list()

for (i in unique(new_predictions$variable)) {
  #i = unique(new_predictions$variable)[1]
  
  df <- filter(new_predictions,
               in_monitoring_area == TRUE,
               variable == i
               )
  
  p[[i]] <- ggplot() +
  geom_sf(data=monitoring_area)  +
  geom_point(data= df, 
            
             aes(x=longitude, y=latitude, col=prediction), alpha=0.3
             ) + 
  facet_wrap(~variable) #+
  #labs(title = "UK-PLS predictions for locations in monitoring area", fill = "")

}

p[1] + theme_bw()

ggarrange(plotlist = p) %>% 
  annotate_figure(top = "UK-PLS predictions for locations in monitoring area")

ggsave(file.path("..", "..", "Kaya", "Block Predictions", "monitoring_predictions.png"),
       height = 11, width = 8)


###########################################################################################
# SAVE PREDICTIONS
###########################################################################################
write.csv(new_predictions, file.path("..", "..", "Kaya", "Block Predictions", "block_predictions.csv"), row.names = F)
saveRDS(new_predictions, file.path("..", "..", "Kaya", "Block Predictions", "block_predictions.rda")) 

###########################################################################################
# TEST - CHECK THAT PREDICTIONS ARE CORRECTLY SAVED
###########################################################################################
test0 <- read.csv(file.path("..", "..", "Kaya", "old", "block_predictions.csv"))
unique(test0$variable)

test <- readRDS(file.path("..", "..", "Kaya", "old", "block_predictions.rda"))
unique(test$variable)
