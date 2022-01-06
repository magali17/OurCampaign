library(data.table)
# --------------------------------------------------------------------------
# Convert string integer to binary - for MA200 and DiscMini
# --------------------------------------------------------------------------
strtobin <- function(number, noBits=8) {
  binary_vector = rev(as.numeric(intToBits(strtoi(number))))
  binary_vector[-(1:(length(binary_vector) - noBits))]
}

# Wrapper that applies strtobin to vector and transposes to get input as rows
binary.flags <- function(data.object, noBits=8){
  t(sapply(data.object, strtobin, noBits=noBits))
}

# --------------------------------------------------------------------------
# NANOSCAN
# --------------------------------------------------------------------------
# Currently assumes subset to ns_status; errors come out of nanoscan as text
# Some messages may have exceeded character limit
# --------------------------------------------------------------------------
ns.status.all <- function(data.object){
  data.object$cpc_flow_err    <- grepl("CPC", data.object$status)
  data.object$charge_flow_err <- grepl("Charge", data.object$status)
  data.object$tilt_err        <- grepl("Tilt", data.object$status)
  data.object$pulse_err       <- grepl("Pulse", data.object$status)
  data.object$temp_err        <- grepl("Temp", data.object$status)
  data.object$inlet_flow_err  <- grepl("Inlet", data.object$status)
  data.object$sheath_flow_err <- grepl("Sheath", data.object$status)
  data.object$reserv_err      <- grepl("Fi", data.object$status)
  data.object$conc_err        <- grepl("Conc", data.object$status)
  data.object$no_err          <- grepl("No errors", data.object$status)
  # Most flow errors are fatal
  data.object$fatal           <- grepl("Sheath", data.object$status) | 
    grepl("CPC", data.object$status) |
    grepl("Inlet", data.object$status)
  # Flag pulse error if > 1 hour consecutive problems
  pulses   <- data.object[grepl("Pulse", data.object$status)]
  uimports <- unique(pulses$import_id)
  for (i in uimports){
    psub     <- pulses[import_id == i]
    np       <- dim(psub)[1]
    if (np > 59){
      tdiff    <- difftime(psub[2:np]$timeint, psub[1:(np-1)]$timeint, units = "mins")
      print(table(as.numeric(tdiff)))
    } }
  data.object[,"value" := NULL]
}

# Column names
ns.flags  <- c("cpc_flow_err","charge_flow_err","tilt_err","pulse_err","temp_err","inlet_flow_err","sheath_flow_err","reserv_err","conc_err")

# --------------------------------------------------------------------------
# DISCMINI
# --------------------------------------------------------------------------
# Discmini codes can be decomposed into binary indicators
# --------------------------------------------------------------------------
disc.status.all <- function(data.object){
  flags  <- rev(c("filter","diffusion","heating","high_temp","high_volt","diffusion_neg","filter_neg","volt_low","dirt","current_low","flow_low","flow_high","filter_zero","diffusion_zero"))
  result <- as.data.frame(binary.flags(data.object[,status], 14))
  setDT(result)
  data.object$no_err <- data.object$status == 0
  data.object[,(flags) := result][,"value" := NULL]
}

# Column names
disc.flags <- rev(c("filter","diffusion","heating","high_temp","high_volt","diffusion_neg","filter_neg","volt_low","dirt","current_low","flow_low","flow_high","filter_zero","diffusion_zero"))

# --------------------------------------------------------------------------
# MA200
# --------------------------------------------------------------------------
# MA200 codes can be decompoesed into binary indicators
# --------------------------------------------------------------------------
ma200.status.all <- function(data.object){
  flags  <- rev(c("none","starting_up","tape_advance","na","optical_sat","timing_err","spot2","unstable_flow","flow_range","manual_time","skip_tape_adv"))
  result <- as.data.frame(binary.flags(data.object[,status], 11))
  setDT(result)
  data.object$no_err <- data.object$status %in% c(0, 512)
  #add new flag columns w/ result as the values
  # why set value to null?? Amanda: "drops the value column, so if it happens to appear in the data object it's removed, so the result can be merged to a table where all the values from the same time stamp are columns (wide instead of long).  I think it's a little harder to reshape the data after converting the flags to binary."
  data.object[,(flags) := result][,"value" := NULL]
}

# Column names
ma200.flags <- rev(c("none","starting_up","tape_advance","na","optical_sat","timing_err","spot2","unstable_flow","flow_range","manual_time","skip_tape_adv"))

# --------------------------------------------------------------------------
# NO2 CAPS
# --------------------------------------------------------------------------
# NO2 status codes have 5 digits and these are categorical
# --------------------------------------------------------------------------
no2.status.all <- function(data.object){
  # 0-pad code to make sure it's 5 digits
  statlen <- sapply(data.object$status, nchar)
  data.object$status[statlen == 1] <- paste0("0000", data.object$status[statlen == 1])
  data.object$status[statlen == 2] <- paste0("000",  data.object$status[statlen == 2])
  data.object$status[statlen == 3] <- paste0("00",   data.object$status[statlen == 3])
  data.object$status[statlen == 4] <- paste0("0",    data.object$status[statlen == 4])
  # Parse to multiple columns
  data.object$pump_on   <- substr(data.object$status, 1, 1) %in% c(1, 3)
  data.object$filter_in <- substr(data.object$status, 1, 1) %in% c(2, 3)
  data.object$bl_none   <- substr(data.object$status, 2, 2) == 0
  data.object$bl_flush  <- substr(data.object$status, 2, 2) == 1
  data.object$bl_meas   <- substr(data.object$status, 2, 2) == 2
  # 1 LED is off (Used only of PMSSA Monitor)
  data.object$led       <- substr(data.object$status, 3, 3) == 1
  # These should all be 0 (NO2 monitor)
  data.object$monitor   <- substr(data.object$status, 4, 4) != 0
  # These should all be 4 - blue (450 nm)
  data.object$wavelength<- substr(data.object$status, 5, 5) != 4
  data.object$no_err    <- data.object$status %in% c("10004","10104")
  data.object[,"value" := NULL]
}

# Column names
no2.flags <- c("pump_on","filter_in","bl_none","bl_flush","bl_meas","led")