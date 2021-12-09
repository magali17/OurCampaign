###########################################################################################
# SETUP
###########################################################################################
# Clear workspace of all objects and unload all extra (non-base) packages
rm(list = ls(all = TRUE))
if (!is.null(sessionInfo()$otherPkgs)) {
  res <- suppressWarnings(
    lapply(paste('package:', names(sessionInfo()$otherPkgs), sep=""),
      detach, character.only=TRUE, unload=TRUE, force=TRUE))
}

pacman::p_load(knitr, kableExtra,
  tidyverse,
  pls, gstat, #variogram()
  sf,
  #VCA, #anovaVCA()
  #broom,
  parallel
  )

set.seed(1)

#functions
source("0_Functions.R")

use_cores <- 4

###########################################################################################
# UPLOAD DATA
###########################################################################################
# mapping variables
lat_long_crs <- 4326  
 
###########################################################################################
# use same train/test set as HEI simulation work. will help w/ reproducibility & consistency
stops0 <- readRDS(file.path("Data", "Output", "stop_data.rda"))

cov_train <- readRDS(file.path("..", "..", "ACT HEI Supp", "act_hei_aim1a", "Output", "mm_cov_train_set.rda"))
cov_test <- readRDS(file.path("..", "..", "ACT HEI Supp", "act_hei_aim1a", "Output", "mm_cov_test_set.rda"))

###########################################################################################
# CALCULATE ANNUAL AVG
###########################################################################################

annual <- stops0 %>%
  group_by(variable, location) %>%
  summarize(
    mean_of_medians = mean(median_value),
    median_of_means = median(mean_value),
    median_of_medians = median(median_value),
    mean_of_means = mean(mean_value)
  ) %>%
  gather("annual", "value", contains(c("mean", "median"))) %>%
  ungroup()

# add covariates
annual_train <- right_join(annual, cov_train)
annual_test <- right_join(annual, cov_test)

#save full dataset for prediction later
keep_names <- names(annual_train)
keep_names <- keep_names[keep_names %in% names(annual_test)]

annual <-rbind(select(annual_train, all_of(keep_names)), select(annual_test, all_of(keep_names)))

saveRDS(annual, file.path("Data", "Output", "annual.rda"))

###########################################################################################
# COMMON VARIABLES
###########################################################################################
cov_names <- cov_train %>%
  select(log_m_to_a1:last_col()) %>%
  names()

pls_comp_n <- 2

#k-folds for CV
k <- 10

#var_names <- unique(annual$variable)

##################################################################################################
# CREATE RANDOM VALIDATION FOLDS FOR TRAINING SET
##################################################################################################
# df = annual_train
# k.=k

add_random_fold <- function(df, k.=k) {
  # #make sure different variables, annual estimates, etc. receive same fold designation
  # set.seed(2)
  
  folds <- df %>%
    distinct(location) %>%
    mutate(random_fold = sample(1:k.,size = nrow(.), replace = T ))
  
  result <- left_join(folds, df)
  
  return(result)
}

############################################################################################################
annual_train <- add_random_fold(annual_train)

############################################################################################################
# UK-PLS FUNCTION 
############################################################################################################
# fn returns UK-PLS predictions. inputs are two spatial objects (simple features). fn automatically transforms these to a lambert projection

# INPUTS
## modeling_data: dataframe that will be used to fit the pls-uk model
## new_data: dataframe with geocovariates for new prediction locations
## cov_names.: name of the modeling covariates
## pls_comp_n.: number of components that should be used to fit PLS models

# OUTPUT
## dataframe with UK predictions at new locations

# modeling_data = group_split(annual_train, variable, annual)[[1]]
# new_data = group_split(annual_test, variable, annual)[[1]]
# cov_names. = cov_names
# pls_comp_n. = pls_comp_n


uk_pls <- function(new_data, modeling_data, cov_names. = cov_names, pls_comp_n. = pls_comp_n) {
  
  lat_long_crs <- 4326
  #lambert projection for UK model
  lambert_proj <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
  ############################################################################################################
  # fit PLS model to estimate fewer components from geocovariates
  pls_model <- plsr(as.formula(paste('value ~', paste(cov_names., collapse = "+"))),
                    data = modeling_data, ncomp = pls_comp_n., scale=T, center=T)
  
  #extract compoent scores for UK
  ## note: can use make.names() to make column names R friendly
  modeling_data_scores <- predict(pls_model, type = "scores") %>% data.frame() %>%  
    # add location & value info
    cbind(modeling_data) %>%
    #cbind(select(modeling_data, variable, location, annual, value, longitude, latitude)) %>%
    #convert back to sf 
    st_as_sf(coords = c('longitude', 'latitude'), crs=lat_long_crs)
  
  new_data_scores <- predict(pls_model, type = "scores", newdata = new_data) %>% data.frame() %>%  
    cbind(new_data) %>%
    #cbind(select(new_data, variable, location, annual, value, longitude, latitude)) %>%
    st_as_sf(coords = c('longitude', 'latitude'), crs=lat_long_crs)
  
  ############################################################################################################
  # fit UK models & predict at new locations
  
  # UK formula 
  uk_formula <- as.formula(paste("value ~", paste0("Comp.", 1:pls_comp_n., collapse = "+")))
  
  # estimate the variogram model: fit a variogram model, offering to the function several different model options (exponential, spherical, and Matern):
  # using lambert coordinates b/c vertical/horizontal units are the same distance
  # the default distance in gstat is 1/3 of the maximum distance (use cutoff option to change this)
  
  v.uk <- variogram(uk_formula, st_transform(modeling_data_scores, lambert_proj) )
  m.uk <- fit.variogram(v.uk, vgm(c("Exp", "Sph", "Mat")) )
  ##make sure Exp/Sph range estimate is at least 0 when little/no correlation in the data 
  m.uk$range[2] <- max(m.uk$range[2], 1)
  #plot(v.uk, m.uk)
  
  # fit UK to the modeling data and predict at the new data locations
  uk_model <- krige(formula = uk_formula, st_transform(modeling_data_scores, lambert_proj), 
                    newdata = st_transform(new_data_scores, lambert_proj), 
                    model = m.uk)
  
  #save predictions
  predictions <- new_data %>% #select(variable, location, annual, value) %>%
    mutate(prediction = uk_model$var1.pred)
  
  return(predictions)
  
}

saveRDS(uk_pls, file.path("Data", "Output", "Predictions", "uk_pls_model.rda"))

##################################################################################################
# CV function
##################################################################################################
# function returns cross-valited predictions for a given dataset. 
# fn works even if folds don't have all numbers in a sequence (e.g., if too few sites)

# INPUTS
## x: a dataframe with which to conduct UK (i.e., a single pollutant, single annual average approach)
# OUTPUTS
## a dataframe with cross-validated predictions

# x = group_split(annual_train, variable, annual)[[1]]
# fold_name = "random_fold"

do_cv <- function (x, fold_name = "random_fold") {
  #code to make sure this fn works even if folds don't have all numbers in a sequence (e.g., if too few sites)
  k = sort(unique(x[[fold_name]]))
  
  df <- data.frame()
  
  for(f in k) {
    #f = k[2]
    modeling_data0 = filter(x, !!as.symbol(fold_name) != f)
    new_data0 = filter(x, !!as.symbol(fold_name) == f)
    
    temp <- uk_pls(new_data = new_data0, modeling_data = modeling_data0) #%>% st_drop_geometry() 
    df <- rbind(df, temp)
  }
  
  return(df)
}

##################################################################################################
# PREDICT
##################################################################################################
common_names <- c("location", "variable", "annual", "value", "prediction", "out_of_sample") 

# cross-validation
cv_predictions <- mclapply(group_split(annual_train, variable, annual), FUN = do_cv, mc.cores = use_cores, 
                         ) %>%
  bind_rows() %>%
  mutate(out_of_sample = "cross-validation") %>%
  select(all_of(common_names))

#test sites
# #UK predictions for the test set, by variable and annual average estimate
#x = group_split(distinct(annual_test, variable, annual), variable, annual)[[1]]

test_predictions <- mclapply(group_split(distinct(annual_test, variable, annual), variable, annual), 
                           function(x) {
                             modeling_data0 = right_join(annual_train, x)
                             new_data0 = right_join(annual_test, x)
                             
                             uk_pls(new_data = new_data0, modeling_data = modeling_data0)
                             }, mc.cores = use_cores
                           ) %>%
  bind_rows() %>%
  mutate(out_of_sample = "test") %>%
  select(all_of(common_names))

# combine predictions
predictions <- rbind(cv_predictions, test_predictions) #%>% select(variable, location, annual, value, prediction, out_of_sample) 

saveRDS(predictions, file.path("Data", "Output", "Predictions", "location_predictions.rda"))

##################################################################################################
# MODEL EVALUATION
##################################################################################################

# dt = group_split(predictions, variable, annual, out_of_sample)[[1]]
# prediction = "prediction"
# reference = "value"

# Fn returns RMSE and MSE-based R2 for a given dataset
validation_stats <- function(dt, prediction = "prediction", reference = "value"){
  
  # MSE of predictions
  MSE_pred <- mean((dt[[reference]] - dt[[prediction]])^2)
  # MSE of observations (for R2 denominator)
  MSE_obs <- mean((dt[[reference]] - mean(dt[[reference]]))^2)
  
  RMSE = sqrt(MSE_pred)
  MSE_based_R2 = max(1 - MSE_pred/MSE_obs, 0)
  reg_based_R2 = cor(dt[[reference]], dt[[prediction]], method = "pearson")^2
  
  result <- distinct(dt, variable, annual, out_of_sample) %>%
    mutate(
      no_sites = nrow(dt),
      RMSE = RMSE,
      MSE_based_R2 = MSE_based_R2,
      reg_based_R2 = reg_based_R2
    )
  
  return(result)
  
}

##################################################################################################

model_performance <- mclapply(group_split(predictions, variable, annual, out_of_sample), validation_stats, mc.cores = 5) %>%
  bind_rows() %>%
  mutate_if(is.numeric, ~round(., 2))


model_performance %>%
  filter(annual == "mean_of_medians") %>%
  select(-c(annual, no_sites)) %>%
  variable_relabel() %>%
  
  kable(caption = "Model performances for annual mean of medians concentrations. N=278 for cross-validation set; N=31 for test set.",
        col.names = c("Pollutant", "Out-of-Sample Set", "RMSE", "MSE-based R2", "Regression-based R2", "UFP Size (nm)")
        ) %>%
  kable_styling()


##################################################################################################
# SAVE FULL DATASET FOR PREDICTION LATER
##################################################################################################

write.csv(model_performance, file.path("Data", "Output", "Predictions", "model_performance.csv"), row.names = F)
 