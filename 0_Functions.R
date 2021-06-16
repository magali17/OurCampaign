

######################################################################################################################
# fn returns colocation plot w/ best fit info, RMSE and MSE-based R2

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
    geom_abline(intercept = 0, slope = 1) +
    #geom_smooth(aes(fill="loess")) +
    geom_smooth(method = "lm", aes(fill="LS")) +
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

add_season <- function(dt, .date_var) {
  
  pacman::p_load(lubridate)
  
  # dt <- aqs_daily
  # .date_var <- "Date.Local"
  
  winter <- "-12-21" #usually winter starts on 21st, sometimes on 22nd 
  spring <- "-03-20"
  summer <- "-06-21" #usually summer starts on 21st, sometimes on 22nd 
  fall <- "-09-23" #usually fall starts on 22nd, sometimes on 23nd. Using 23rd for 2019 mobile monitoring campaign 
  
  dt <- dt %>%
    rename(date_var = .date_var) %>%
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
# * function for a standard boxplot with different whisker definitions to avoid plotting extreme/outlier points
#function takes df and returns summary statistics for plotting alternative boxplots with quantiles: 10, 25, 50, 75 and 90. this reduces the plotting of outliers, which are typically problematic when dealign with large datasets. 

alt_boxplot <- function(df, var, min_q=0.025, max_q=0.975){
  df <- df %>%
    rename(var = var) %>%
    
    #calculate quantiles
    summarize(
      Qmin = quantile(var, min_q),
      Q25 = quantile(var, 0.25),
      Q50 = quantile(var, 0.50),
      Q75 = quantile(var, 0.75),
      Qmax = quantile(var, max_q)
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


######################################################################################################################

# fn relabels variable (pollutants) to include units and a clearer name for presentation purposes

variable_relabel <- function(dt, var = "variable") {
  dt <- dt %>%
    rename(var=var) %>%
    mutate(
      ufp_range_nm = ifelse(var=="ns_total_conc", "10-420",
                         ifelse(var=="pmdisc_number", "10-700",
                                ifelse(var=="pnc_noscreen", "20-1,000",
                                       ifelse(var=="pnc_screen", "50-1,000", "-"
                                       )))),
      #ufp_range_nm = factor(ufp_range_nm, levels = c("10-700", "10-420", "20-1,000", "50-1,000")),
      
      var = recode_factor(factor(var),
                          "co_ppm" = "CO (ppm)",
                          "co2_umol_mol" = "CO2 (ppm)",
                          "ma200_ir_bc1" = "BC (ng/m3)",
                          "neph_bscat" = "Neph (bscat/m)",
                          "pm2.5_ug_m3" = "PM2.5 (ug/m3)",
                          "no2" = "NO2 (ppb)",
                          "ns_total_conc" = "UFP (pt/cm3)", #, 10-420 nm (pt/cm3)", #"UFP_NanoScan (pt/cm3)",
                          "pmdisc_number" = "UFP (pt/cm3)", #, 10-700 nm (pt/cm3)", #"UFP_DiSCmini (pt/cm3)",
                          "pnc_noscreen" = "UFP (pt/cm3)", #, 20-1,000 nm (pt/cm3)", #"UFP_PTRAK (pt/cm3)",
                          "pnc_screen" = "UFP (pt/cm3)", #, 50-1,000 nm (pt/cm3)" #"UFP_PTRAK_screened (pt/cm3)"
      ),
      var = factor(var, levels = c("Neph (bscat/m)",  "PM2.5 (ug/m3)",
                                   "BC (ng/m3)", 
                                   "UFP (pt/cm3)", #, 10-700 nm (pt/cm3)", "UFP, 10-420 nm (pt/cm3)", "UFP, 20-1,000 nm (pt/cm3)", "UFP, 50-1,000 nm (pt/cm3)",
                                   "NO2 (ppb)",
                                   "CO2 (ppm)", 
                                   
                                   "CO (ppm)"
                                   )
                   
                   ),
      
      # #drop everything after a space or comma
      # pollutant = gsub(" .*|,.*", "", var),
      # pollutant = factor(pollutant, levels = c("")
      #                      )
    )
  
  names(dt)[names(dt) == "var"] <- var
  
  return(dt)
  
}

######################################################################################################################

