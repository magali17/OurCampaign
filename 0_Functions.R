# geocovariate preprocessing
######################################################################################################################
# function labels covariates 

split_cov_name <- function(dt, cov) {
  dt <- suppressWarnings(dt %>%
                           rename(cov_full_name = cov) %>%
                           mutate(
                             Buffer = substr(cov_full_name, nchar(cov_full_name)-4, nchar(cov_full_name)),
                             #for non-buffered covariates, use "0"
                             Buffer = as.numeric(ifelse(!is.na(as.integer(Buffer)), Buffer, 0)),
                             #drop buffer repetition
                             cov = ifelse(Buffer==0, cov_full_name, substr(cov_full_name, 1, nchar(cov_full_name)-5) )
                           )
  ) %>%
    select(contains("cov"), Buffer, everything())
  
  # elevation
  dt$cov[grepl("^elev_.+_above$", dt$cov_full_name)] <- "elev_above"
  dt$cov[grepl("^elev_.+_below$", dt$cov_full_name)] <- "elev_below"
  dt$cov[grepl("^elev_.+_stdev$", dt$cov_full_name)] <- "elev_stdev"
  dt$cov[grepl("^elev_.+_at_elev$", dt$cov_full_name)] <- "elev_at_elev"
  
  dt$Buffer[grepl("_1k_", dt$cov_full_name)] <- 1000
  dt$Buffer[grepl("_5k_", dt$cov_full_name)] <- 5000
  
  return(dt)
}
######################################################################################################################

#cov_name = "covariate"

covariate_table_fn <- function(dt, cov_name) {
  
  dt[[cov_name]] = as.character(dt[[cov_name]])
  
  dt %>% 
    #mutate(covariate= as.character(!!as.symbol(cov_name) )) %>%
    split_cov_name(cov = cov_name)%>%
    group_by(Covariate = cov) %>%
    dplyr::summarize("Buffers" = paste0(sort(unique(Buffer)), collapse = ", ")) %>%
    mutate(Description = recode(Covariate,
                                bus_s = "sum of bus routes",
                                log_m_to_bus = "log meters to bus route",
                              elev_above = "number of points (out of 24) more than 20 m and 50 m uphill of a location for a 1000 m and 5000 m buffer, respectively",
                              elev_below = "number of points (out of 24) more than 20 m and 50 m downhill of a location for a 1000 m and 5000 m buffer, respectively",
                              elev_at_elev = "number of points (out of 24) within 20 m and 50 m of the location' elevation for a 1000 m and 5000 m buffer, respectively",
                              elev_stdev = "standard deviation of elevation of 20 points surrounding the location",
                              elev_elevation = "elevation above sea level in meters",
                              imp_a = "average imperviousness",
                              intersect_a1_a3_s = "number of A1-A3 road intersections",
                              intersect_a2_a2_s = "number of A2-A2 road intersections",
                              intersect_a2_a3_s = "number of A2-A3 road intersections",
                              intersect_a3_a3_s = "number of A3-A3 road intersections",
                              ll_a1_s = "length of A1 roads",
                              ll_a2_s = "length of A2 roads",
                              ll_a3_s = "length of A3 roads",
                              ll_a23_s = "length of A2 and A3 roads",
                              log_m_to_a1 = "log meters to closest A1 road",
                              log_m_to_a1_a1_intersect = "log meters to closest A1-A1 road intersection",
                              log_m_to_a1_a2_intersect= "log meters to closest A1-A2 road intersection",
                              log_m_to_a1_a3_intersect = "log meters to closest A1-A3 road intersection",
                              log_m_to_a123 = "log meters to closest A1, A2 or A3 road",
                              log_m_to_a2 = "log meters to closest A2 road",
                              log_m_to_a2_a2_intersect = "log meters to closest A2-A2 road intersection",
                              log_m_to_a2_a3_intersect= "log meters to closest A2-A3 road intersection",
                              log_m_to_a23 = "log meters to closest A2 or A3 road",
                              log_m_to_a3 = "log meters to closest A3 road",
                              log_m_to_a3_a3_intersect = "log meters to closest A3-A3 road intersection",
                              log_m_to_airp = "log meters to closest airport",
                              log_m_to_coast = "log meters to closest coastline",
                              log_m_to_comm = "log meters to closest commercial and services area",
                              log_m_to_l_airp = "log meters to closest large airport",
                              log_m_to_l_port = "log meters to closest large port",
                              log_m_to_m_port = "log meters to closest medium port",
                              log_m_to_rr = "log meters to closest railroad",
                              log_m_to_ry = "log meters to closest rail yard",
                              log_m_to_truck = "log meters to closest truck route",
                              log_m_to_waterway = "log meters to closest waterway",
                              ndvi_q25_a = "NDVI (25th quantile)",
                              ndvi_q50_a = "NDVI (50th quantile)",
                              ndvi_q75_a = "NDVI (75th quantile)",
                              ndvi_summer_a = "average summer time NDVI",
                              ndvi_winter_a = "average winter time NDVI",
                              pop_s = "2000 population density",
                              pop10_s = "2010 population density",
                              pop90_s = "1990 population density",
                              rlu_barren_p = "proportion of barren land",
                              rlu_decid_forest_p = "proportion of deciduous forest",
                              rlu_dev_hi_p = "proportion of highly developed land (e.g., commercial and services; industrial; transportation, communication and utilities)",
                              rlu_dev_lo_p = "proportion of low developed land  (e.g., residential)",
                              rlu_dev_med_p = "proportion of medium developed land (e.g., residential)",
                              rlu_dev_open_p = "proportion of developed open land",
                              rlu_evergreen_p = "proportion of evergreen forest",
                              rlu_herb_wetland_p = "proportion of herb (nonforested) wetland",
                              rlu_mix_forest_p = "proportion of mixed forest",
                              rlu_water_p = "proportion of water",
                              rlu_woody_wetland_p = "proportion of woody wetland",
                              tl_s = "length of truck routes",
                              
                              emissions_s = "sum of emissions (g) on highways and major roads",
                              
                              em_CO_s = "sum of CO stack emissions",
                              em_NOx_s = "sum of NOx stack emissions",
                              em_PM10_s = "sum of PM10 stack emissions",
                              em_PM25_s = "sum of PM2.5 stack emissions",
                              em_SO2_s = "sum of SO2 stack emissions",
                              no2_behr = "columnar NO2, mean from 2005-2007",
                              no2_behr_2005 = "columnar NO2 in 2005",
                              no2_behr_2006 = "columnar NO2 in 2006",
                              no2_behr_2007 = "columnar NO2 in 2007",
                              
  ),
  Kind = ifelse(grepl("elev", Covariate), "elevation", 
                ifelse(grepl("^imp", Covariate), "imperviousness",
                       ifelse(grepl("_a[1-3]", Covariate), "roads",
                              ifelse(grepl("airp", Covariate), "airports",
                                     ifelse(grepl("coast", Covariate), "coast",
                                            ifelse(grepl("comm", Covariate), "commercial and services",
                                                   ifelse(grepl("port", Covariate), "port",
                                                          ifelse(grepl("rr$|ry$", Covariate), "railroads, rail yards",
                                                                 ifelse(grepl("truck|tl_", Covariate), "truck routes",
                                                                        # --> ? rlu also has a water variable. keep as "land use"?
                                                                        ifelse(grepl("water", Covariate), "water",
                                                                               ifelse(grepl("ndvi", Covariate), "NDVI",
                                                                                      ifelse(grepl("^pop", Covariate), "population",
                                                                                             ifelse(grepl("^rlu", Covariate), "land use",
                                                                                                    
                                                                                                    ### --> UPDATE emissions text
                                                                                                    ifelse(grepl("emissions|^em_", Covariate), "stack emissions", 
                                                                                                           ifelse(grepl("^no2_behr", Covariate), "columnar NO2",
                                                                                                                  ifelse(grepl("bus", Covariate), "bus", "other")
                                                                                                    ) ) ) ) )))))))))))
  ) %>%
  select(Kind, everything()) %>%
  arrange(Kind)


}

######################################################################################################################
# fn returns colocation plot w/ best fit info, RMSE and MSE-based R2
## blue line is linear fit

colo.plot <- function(data.wide=mm.wide, 
                      x.variable, x.label = "",
                      y.variable, y.label = "",
                      col.by = "",
                      shape.var="",
                      alpha_value = 0.3,
                      mytitle = "", title_width = 60,
                      mysubtitle = NULL,
                      mycaption = NULL,
                      int_digits = 0,
                      slope_digits = 2,
                      r2.digits = 2, 
                      rmse.digits = 0,
                      convert_rmse_r2_to_native_scale=F
) {
  
  #if label is left blank, use variable name
  if(x.label == ""){x.label <- x.variable}
  if(y.label == ""){y.label <- y.variable}
  
  data.wide <- data.wide %>% drop_na(x.variable, y.variable)  
  
  lm1 <- lm(formula(paste(y.variable, "~", x.variable)), data = data.wide)
  #summary(lm1)
  
  
  ################################################
  ## ?? need this fns inside this fn???
  #returns MSE
  mse <- function(obs, pred){
    mean((obs - pred)^2)
  }
  
  rmse <- function(obs, pred){
    sqrt(mean((obs - pred)^2))
  }
  
  #returns MSE-based R2
  r2_mse_based <- function(obs, pred) {
    mse.est <- mse(obs, pred)
    r2 <- 1- mse.est/mean((obs - mean(obs))^2)
    max(0, r2)
  }  
  
  ################################################ 
  
  
  #rmse
  if (convert_rmse_r2_to_native_scale==T) {
    rmse <- rmse(obs = exp(data.wide[[x.variable]]), pred = exp(data.wide[[y.variable]])) %>% 
      round(digits = rmse.digits)
    
    r2 <- r2_mse_based(obs = exp(data.wide[[x.variable]]), pred = exp(data.wide[[y.variable]])) %>%
      round(r2.digits)
  } 
  
  else {
    rmse <- rmse(obs = data.wide[[x.variable]], pred = data.wide[[y.variable]]) %>% 
      round(digits = rmse.digits)
    
    r2 <- r2_mse_based(obs = data.wide[[x.variable]], pred = data.wide[[y.variable]]) %>%
      round(r2.digits)
  }
  
  
  fit.info <- paste0("y = ", round(coef(lm1)[1], int_digits), " + ", round(coef(lm1)[2], slope_digits), 
                     "x \nMSE-R2 = ", r2,  
                     "\nRMSE = ", rmse,
                     "\nNo. Pairs = ", nrow(data.wide))
  
  max_plot <- max(max(data.wide[[x.variable]]), max(data.wide[[y.variable]]) )
  min_plot <- min(min(data.wide[[x.variable]]), min(data.wide[[y.variable]]) )
  
  #compare  
  p <- data.wide %>%
    ggplot(aes(x= data.wide[[x.variable]], y= data.wide[[y.variable]])) + 
    geom_point(alpha=alpha_value, aes(col = data.wide[[col.by]],
                                      shape = data.wide[[shape.var]]
    )) + 
    coord_fixed() +
    geom_abline(intercept = 0, slope = 1, linetype=2, alpha=0.5) +
    #geom_smooth(aes(fill="loess")) +
    geom_smooth(method = "lm", se = F) +
    xlim(min_plot, max_plot) +
    ylim(min_plot, max_plot) +
    labs(title = #wrapper(
           mytitle, width = title_width
         #)
         ,
         subtitle = mysubtitle,
         caption = mycaption,
         x = x.label,
         y = y.label,
         col = col.by, 
         shape = shape.var,
         fill = "fit") +
    annotate("text", -Inf, Inf, label = fit.info, hjust = 0, vjust = 1)
  
  return(p)
  
}


######################################################################################################################

# fn adds sesason to a given dataset with a date variable. Uses typical equinox/solstice dates
#dt = test
#.date_var = "date"


add_season <- function(dt, .date_var) {
  
  pacman::p_load(lubridate)
  
  # dt <- aqs_daily
  # .date_var <- "Date.Local"
  
  winter <- "-12-21" #usually winter starts on 21st, sometimes on 22nd 
  spring <- "-03-20"
  summer <- "-06-21" #usually summer starts on 21st, sometimes on 22nd 
  fall <- "-09-23" #usually fall starts on 22nd, sometimes on 23nd. Using 23rd for 2019 mobile monitoring campaign 
  
  dt <- dt %>%
    dplyr::rename(date_var = .date_var) %>%
    #make sure variable is in date format
    mutate(date_var = as.Date(date_var),
           season = factor(ifelse((date_var >= ymd(paste0((year(date_var)-1), winter)) & date_var < ymd(paste0(year(date_var), spring))) |
                                    date_var >= ymd(paste0(year(date_var), winter)), "winter",
                                  ifelse(date_var >= ymd(paste0(year(date_var), spring)) &
                                           date_var < ymd(paste0(year(date_var), summer)), "spring",
                                         ifelse(date_var >= ymd(paste0(year(date_var), summer)) &
                                                  date_var < ymd(paste0(year(date_var), fall)), "summer", 
                                                ifelse( date_var >= ymd(paste0(year(date_var), fall)) &
                                                          date_var < ymd(paste0(year(date_var), winter)), "fall", 
                                                        NA)))), 
                           levels = c("spring", "summer", "fall", "winter"))
    )
  
  #change time variable back to what it was originally
  names(dt)[names(dt) %in% "date_var"] <- .date_var
  
  return(dt)
  
}

######################################################################################################################
# windsorizes a value
# value = "median_value"
# dt = stops0

winsorize_fn <- function(dt, value, trim_quantile =0.05) {

  #dt <- rename(dt, value = value)
  
  dt1 <- dt %>%
  group_by(variable, location) %>%
    mutate(
    # wind_median_value = ifelse(median_value == max(median_value), max(median_value[median_value!=max(median_value)]),
    #                            ifelse(median_value == min(median_value), min(median_value[median_value!=min(median_value)]),
    #                                   median_value
    #                            )),
    # win_value = ifelse(value > quantile(value, 1-trim_quantile), quantile(value, 1-trim_quantile), 
    #                              ifelse(value < quantile(value, trim_quantile), quantile(value, trim_quantile),
    #                                     value))
      win_value = ifelse(!!as.symbol(value) > quantile(!!as.symbol(value), 1-trim_quantile), quantile(!!as.symbol(value), 1-trim_quantile), 
                         ifelse(!!as.symbol(value) < quantile(!!as.symbol(value), trim_quantile), quantile(!!as.symbol(value), trim_quantile),
                                !!as.symbol(value)))
    )
  
  #names(dt1)[names(dt1) == value] <- value
  
  
  }


######################################################################################################################
# * function for a standard boxplot with different whisker definitions to avoid plotting extreme/outlier points
#function takes df and returns summary statistics for plotting alternative boxplots with quantiles: 10, 25, 50, 75 and 90. this reduces the plotting of outliers, which are typically problematic when dealign with large datasets. 

alt_boxplot <- function(df, var, min_q=0.025, max_q=0.975){
  df <- df %>%
    rename(var = var) %>%
    
    #calculate quantiles
    summarize(
      N = n(),
      Min = min(var),
      Qmin = quantile(var, min_q),
      Q25 = quantile(var, 0.25),
      Q50 = quantile(var, 0.50),
      Q75 = quantile(var, 0.75),
      Qmax = quantile(var, max_q),
      Max = max(var)
    )
  
  names(df)[names(df)==var] <- var
  
  return(df) 
  
}


######################################################################################################################

#fn returns MSE
mse <- function(obs, pred){
  mean((obs - pred)^2)
}

rmse <- function(obs, pred){
  sqrt(mean((obs - pred)^2))
}


#fn returns MSE-based R2
r2_mse_based <- function(obs, pred) {
  mse.est <- mse(obs, pred)
  r2 <- 1- mse.est/mean((obs - mean(obs))^2)
  max(0, r2)
}  


######################################################################################################################

# facet_wrap_equal() function acts like facet_wrap() in ggplot but it sets the axes ranges (min/max) of each facet to same scale so that the 1-1 line is always down the middle :D !!

# code source: https://fishandwhistle.net/post/2018/modifying-facet-scales-in-ggplot2/ 

FacetEqualWrap <- ggproto(
  "FacetEqualWrap", FacetWrap,
  
  train_scales = function(self, x_scales, y_scales, layout, data, params) {
    
    # doesn't make sense if there is not an x *and* y scale
    if (is.null(x_scales) || is.null(x_scales)) {
      stop("X and Y scales required for facet_equal_wrap")
    }
    
    # regular training of scales
    ggproto_parent(FacetWrap, self)$train_scales(x_scales, y_scales, layout, data, params)
    
    # switched training of scales (x and y and y on x)
    for (layer_data in data) {
      match_id <- match(layer_data$PANEL, layout$PANEL)
      
      x_vars <- intersect(x_scales[[1]]$aesthetics, names(layer_data))
      y_vars <- intersect(y_scales[[1]]$aesthetics, names(layer_data))
      
      SCALE_X <- layout$SCALE_X[match_id]
      ggplot2:::scale_apply(layer_data, y_vars, "train", SCALE_X, x_scales)
      
      SCALE_Y <- layout$SCALE_Y[match_id]
      ggplot2:::scale_apply(layer_data, x_vars, "train", SCALE_Y, y_scales)
    }
    
  }
)

facet_wrap_equal <- function(...) {
  # take advantage of the sanitizing that happens in facet_wrap
  facet_super <- facet_wrap(...)
  
  ggproto(NULL, FacetEqualWrap,
          shrink = facet_super$shrink,
          params = facet_super$params
  )
}


#same as above but for facet_grid()
FacetEqualGrid <- ggproto(
  "FacetEqualGrid", FacetGrid,
  
  train_scales = function(self, x_scales, y_scales, layout, data, params) {
    
    # doesn't make sense if there is not an x *and* y scale
    if (is.null(x_scales) || is.null(x_scales)) {
      stop("X and Y scales required for facet_equal_wrap")
    }
    
    # regular training of scales
    ggproto_parent(FacetGrid, self)$train_scales(x_scales, y_scales, layout, data, params)
    
    # switched training of scales (x and y and y on x)
    for (layer_data in data) {
      match_id <- match(layer_data$PANEL, layout$PANEL)
      
      x_vars <- intersect(x_scales[[1]]$aesthetics, names(layer_data))
      y_vars <- intersect(y_scales[[1]]$aesthetics, names(layer_data))
      
      SCALE_X <- layout$SCALE_X[match_id]
      ggplot2:::scale_apply(layer_data, y_vars, "train", SCALE_X, x_scales)
      
      SCALE_Y <- layout$SCALE_Y[match_id]
      ggplot2:::scale_apply(layer_data, x_vars, "train", SCALE_Y, y_scales)
    }
    
  }
)

facet_grid_equal <- function(...) {
  # take advantage of the sanitizing that happens in facet_wrap
  facet_super <- facet_grid(...)
  
  ggproto(NULL, FacetEqualGrid,
          shrink = facet_super$shrink,
          params = facet_super$params
  )
}

######################################################################################################################

# fn relabels variable (pollutants) to include units and a clearer name for presentation purposes

# dt = grid_predictions2 
# var = "variable"
# reverse_labels = TRUE

variable_relabel <- function(dt, var = "variable", keep_original_var = FALSE, 
                             # if original labels don't exist, add them
                             reverse_labels = FALSE) {
  
  dt <- dt %>% 
    dplyr::rename(var=var) 
  
  if(reverse_labels == FALSE) {
    
    if(keep_original_var==TRUE){dt$var0 <- dt$var}
    
    dt <- dt %>%
      mutate(
        ufp_range_nm = ifelse(var=="ns_total_conc", "10-420 nm",
                           ifelse(var=="pmdisc_number", "10-700 nm",
                                  ifelse(var=="pnc_noscreen", "20-1,000 nm",
                                         ifelse(var=="pnc_screen", "36-1,000 nm", 
                                                ifelse(var=="pnc_20_36", "20-36 nm",
                                                       "-"
                                         ))))),
        ufp_instrument = ifelse(var=="ns_total_conc", "NanoScan",
                              ifelse(var=="pmdisc_number", "DiSCmini",
                                     ifelse(var=="pnc_noscreen", "P-TRAK",
                                            ifelse(var=="pnc_screen", "Screened P-TRAK",
                                                   "-"
                                                   )))),
        ufp_instrument = factor(ufp_instrument, levels = c("P-TRAK", "Screened P-TRAK", "NanoScan", "DiSCmini",  "-", "NA")),
        
        var = recode_factor(factor(var),
                            "co_ppm" = "CO (ppm)",
                            "co2_umol_mol" = "CO2 (ppm)",
                            "ma200_ir_bc1" = "BC (ng/m3)",
                            "neph_bscat" = "Neph (bscat/m)",
                            "pm2.5_ug_m3" = "PM2.5 (ug/m3)",
                            "no2" = "NO2 (ppb)",
                            "pnc_noscreen" = "PNC (pt/cm3)", #, 20-1,000 nm
                            "pnc_screen" = "PNC (pt/cm3)", #, 50-1,000 nm  
                            "pnc_20_36" = "PNC (pt/cm3)", #, 20-50 nm  
                            "ns_total_conc" = "PNC (pt/cm3)", #, 10-420 nm  
                            "pmdisc_number" = "PNC (pt/cm3)", #, 10-700 nm  
                            
        ),
        var = factor(var, levels = c("PNC (pt/cm3)", 
                                     "BC (ng/m3)", 
                                     "NO2 (ppb)",
                                     "Neph (bscat/m)",  "PM2.5 (ug/m3)",
                                     "CO2 (ppm)", 
                                     "CO (ppm)"
                                     )
                     )
      )
    
    }
  
  # reverse labels
  # make labels R friendly 
  if(reverse_labels == TRUE) {
    
    dt <- mutate(dt,
      ufp_range_nm = as.factor(ufp_range_nm),
      
      variable0 = case_when(
      var == "CO (ppm)" ~ "co_ppm",
      var == "CO2 (ppm)" ~ "co2_umol_mol",
      var == "BC (ng/m3)" ~ "ma200_ir_bc1",
      var == "PM2.5 (ug/m3)" ~ "pm2.5_ug_m3",
      var == "NO2 (ppb)" ~ "no2",
      ufp_range_nm == "20-1,000 nm" ~ "pnc_noscreen",
      ufp_range_nm == "36-1,000 nm" ~ "pnc_screen",
      ufp_range_nm == "10-420 nm" ~ "ns_total_conc",
      ufp_range_nm == "10-700 nm" ~ "pmdisc_number",
      TRUE ~ ""
      )
      )

  }
  
  
  names(dt)[names(dt) == "var"] <- var
  
  return(dt)
  
} 

######################################################################################################################

