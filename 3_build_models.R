###########################################################################################
# NOTES
###########################################################################################
# mclapply is used in this code, but it this is not supported on windows when mc.cores >1. 
# windows users may want to consider replacing mclapply() with lapply() instead

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

image_path <- file.path("..", "Manuscript", "Images")

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
  #winsorize median values 
  winsorize_fn(value = "median_value") %>%  
  
  summarize(
    mean_of_medians = mean(median_value),
    median_of_means = median(mean_value),
    median_of_medians = median(median_value),
    mean_of_means = mean(mean_value),
    
    #mean_of_wind_medians = mean(wind_median_value),
    mean_of_win_medians = mean(win_value),
    #mean_of_q05_win_medians = mean(q05win_median_value)
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

# for sharing 
annual %>%
  select(variable, location, longitude, latitude, annual, value) %>%
  write.csv(., file.path("Data", "Output", "annual.csv"), row.names = F)


###########################################################################################
# TEST - see influential points
###########################################################################################

# high_annual_sites <- annual %>%
#   filter(annual == "mean_of_medians",
#          grepl("pnc|pmdisc|ns_total", variable)
#          ) %>%
#   group_by(variable)%>%
#   filter(value > quantile(value, 0.98)) %>% 
#   arrange(variable, desc(value)) %>%
#   rename(annual_value = value)
#   #select(variable, location, annual_value = value)
# 
# #MS0229
# 
# # stops0 %>%
# #   filter(grepl("pmdisc", variable),
# #          #grepl("pnc|pmdisc|ns_total", variable),
# #          location %in% c("MS0229")
# #          ) %>% View()
# 
# 
# # most sites w/ high UFP readings are influenced by 1+ very high reading
# high_annual_sites %>%
#   left_join(stops0) %>% # View()
# 
#   ggplot(aes(y=location, x=median_value, col=annual_value)) +
#   facet_wrap(~variable, scales="free") +
#   geom_point(alpha=0.3) +
#   labs(x = "2-min Median",
#        col = "Annual mean\n of median",
#        title = "2-min medians for sites with high annual averages"
#        )
# ggsave(file.path(image_path, "Temp", "Windsorized", "stop_medians_of_high_annual_avgs.png"), height = 6, width = 6)


###########################################################################################
# COMMON VARIABLES
###########################################################################################
cov_names <- cov_train %>%
  select(log_m_to_a1:last_col()) %>%
  names()

pls_comp_n <- 3 #note, that PM2.5 & no2 did a tiny bit better w/ 4 PLS components
save(pls_comp_n, file = file.path("Data", "Output", "Objects", "pls_comp_n.rda"))

#k-folds for CV
k <- 10

##################################################################################################
# CREATE RANDOM VALIDATION FOLDS FOR TRAINING SET
##################################################################################################

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


uk_pls <- function(new_data, modeling_data, cov_names. = cov_names, pls_comp_n. = pls_comp_n, fn_result = "predictions") {
  
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
  
  # fit UK to the modeling data and predict at the new data locations
  uk_model <- krige(formula = uk_formula, st_transform(modeling_data_scores, lambert_proj), 
                    newdata = st_transform(new_data_scores, lambert_proj), 
                    model = m.uk)
  
  #save predictions
  predictions <- new_data %>% #select(variable, location, annual, value) %>%
    mutate(prediction = uk_model$var1.pred)
  
  if(fn_result == "predictions") {return(predictions)}
  if(fn_result == "models") {
    result = list(
      variable = first(modeling_data$variable),
      annual = first(modeling_data$annual),
      pls_model = pls_model, 
      variogram_model = m.uk
      )
    return(result)
    }
  
  
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

do_cv <- function (x, fold_name = "random_fold") {
  #code to make sure this fn works even if folds don't have all numbers in a sequence (e.g., if too few sites)
  k = sort(unique(x[[fold_name]]))
  
  df <- data.frame()
  
  for(f in k) {
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

#same as above but returns the pls & variogram model used in UK with the training-validation set
train_validation_models <- mclapply(group_split(distinct(filter(annual_test, annual=="mean_of_medians"), variable, annual), variable, annual), 
                             function(x) {
                               modeling_data0 = right_join(annual_train, x)
                               new_data0 = right_join(annual_test, x)
                               
                               uk_pls(new_data = new_data0, modeling_data = modeling_data0, fn_result = "models")
                             }, mc.cores = use_cores
                             )

saveRDS(train_validation_models, file.path("Data", "Output", "Predictions", "train_validation_models.rda"))


# combine predictions
predictions <- rbind(cv_predictions, test_predictions) 

saveRDS(predictions, file.path("Data", "Output", "Predictions", "monitoring_location_predictions.rda"))

## for sharing
predictions %>%
  left_join(select(annual, location, longitude, latitude)) %>% 
  select (location, longitude, latitude, variable, annual, value, prediction) %>%# View()
  write.csv(., file.path("Data", "Output", "Predictions", "monitoring_location_estimates_and_predictions.csv"), row.names = F)


##################################################################################################
# MODEL EVALUATION
##################################################################################################

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

# table
model_performance %>%
  filter(annual %in% c("mean_of_medians", "mean_of_win_medians")) %>%
  select(-c(no_sites)) %>%
  variable_relabel() %>%
  select(-ufp_range_nm) %>%
  select(c(variable, ufp_instrument, everything())) %>%

  kable(caption = "Model performances for annual mean of medians concentrations. N=278 for cross-validation set; N=31 for test set.",
        col.names = c("Pollutant", "PNC Instrument", "Annual", "Out-of-Sample Set", "RMSE", "MSE-based R2", "Regression-based R2")
        ) %>%
  kable_styling()


##################################################################################################
# TEST - COMPARE MODEL PERFORMANCES

# model_performance %>%
#   filter(annual %in% c("mean_of_medians", "mean_of_win_medians")) %>%
#   variable_relabel() %>%
#   gather(performance, value, MSE_based_R2, RMSE) %>%
#   
#   ggplot(aes(x=out_of_sample, y=value, col=annual, shape=ufp_range_nm)) + 
#     facet_wrap(~variable+performance, scales="free") + 
#     geom_point() + 
#   labs(title = "model performance of mean of medians vs mean of windsorized medians")
# 
# ggsave(file.path(image_path, "Temp", "Windsorized", "model_performance.png"), height = 8, width = 11)
  
  
##################################################################################################
# SAVE FULL DATASET FOR PREDICTION LATER
##################################################################################################

model_performance %>%
  write.csv(., file.path("Data", "Output", "Predictions", "model_performance.csv"), row.names = F)




 

 