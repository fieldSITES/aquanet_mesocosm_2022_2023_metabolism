# Fit LakeMetabolizer models to Mesocosm experiment data

rm(list = ls())
graphics.off()

Sys.setenv(TZ = "UTC")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(data.table)
library(LakeMetabolizer)
library(stringr)
library(lubridate)
library(ggplot2)

source("compile_wind_ts.R")
source("compile_sensor_ts.R")

##### Settings -----
folder_data = "../SensorData/"
folder_wind = "../WindData/"
folder_out = "../Output/"
folder_plots = file.path(folder_out, "Plots")

remove_days_with_na = F # If T, model is fitted on all days. If not, NAs are returned on days with NAs
gap_filling = T # If T, some minor gap-filling is applied. Larger gaps are left.
max_gap = 40 * 60 # Linear interpolation if gaps are up to this duration, provided in seconds
suffix = ""

years = c(2022, 2023)
sites = c("Erken", "Bolmen", "Skogaryd")
lat_sites = c(Erken = 59.8, Bolmen = 56.9, Skogaryd = 58.4)
elev_sites = c(Erken = 12.0, Bolmen = 142.0, Skogaryd = 74) # Info from https://minkarta.lantmateriet.se/

wind_reduction = 0.5 # What fraction of wind to use (due to sheltering)

depth_meso = 1.0 # Depth of mesocosm (m)

function_to_fit = metab.mle # metab.mle or metab.kalman

# Quality check
apply_gpp_r_bounds = T # Use min_gpp and max_r
min_gpp = -0.5 # Values lower than this are set to NA
max_r = 0.5 # Values higher than this are set to NA
sensors_to_remove = list(Bolmen = list("2023" = "01")) # Number at the end refers to sensor number

fig_type = "png"
fig_width = 8
fig_height = 6

##### Fit metabolism models -----
if(!dir.exists(folder_plots)){
  dir.create(folder_plots, recursive = T)
}

# Progress bar
progressBar = txtProgressBar(min = 0, max = length(years) * length(sites), initial = 0)
progress = 0

for(i in years){
  for(j in sites){
    if(i == 2022 & j == "Skogaryd"){
      progress = progress + 1
      setTxtProgressBar(progressBar, progress)
      next
    }
    
    # Read meteorology - Reduce wind by certain fraction to account of sheltering
    # effect of mesocosms. 
    df_met = compile_wind_ts(folder_wind, site = j, year = i)
    df_met[, wnd := wnd * wind_reduction]
    
    # Read sensor data
    df_exp = compile_sensor_ts(folder_data, site = j, year = i)
    
    meas_frequency_time_step = difftime(unique(df_exp$datetime)[2],
                                        unique(df_exp$datetime)[1],
                                        units = "secs")
    
    # Interpolate wind speed to get same resolution as sensor data
    df_met = merge(df_met, list(datetime = seq.POSIXt(df_met[1L, datetime],
                                                      df_met[.N, datetime],
                                                      by = meas_frequency_time_step)),
                   by = "datetime", all.y = T)
    df_met[, wnd := approx(datetime, wnd, xout = datetime)$y]
    
    # Gap-filling of sensor data
    if(gap_filling){
      max_gap_num = ceiling(max_gap / as.numeric(meas_frequency_time_step))
      
      meas_cols = names(df_exp)[-c(1, 2, 3)]
      
      for(sens in unique(df_exp$Sensor)){
        df_temp = df_exp[Sensor == sens]
        
        for(the_col in meas_cols){
          df_rle = rle(!is.na(df_temp[, ..the_col][[1]]))
          setDT(df_rle)
          df_rle[, cumlength := cumsum(lengths)]
          
          for(k in seq_len(nrow(df_rle))[-1]){
            if(df_rle[k, values == FALSE & df_rle[k, lengths] <= max_gap_num]){
              ind_start = df_rle[k - 1, cumlength]
              ind_end = df_rle[k, cumlength] + 1
              if(ind_end > nrow(df_temp)) next
              
              df_temp[ind_start:ind_end, (the_col) := approx(datetime, get(the_col), xout = datetime)$y]
            }
          }
        }
        
        # After doing all columns, merge back into df_exp
        df_exp[Sensor == sens, (meas_cols) := df_temp[, meas_cols, with = FALSE]]
      }
      
      # Now do the same for wind speed
      df_rle2 = rle(!is.na(df_met[, wnd]))
      setDT(df_rle2)
      df_rle2[, cumlength := cumsum(lengths)]
      for(k in seq_len(nrow(df_rle2))[-1]){
        if(df_rle2[k, values == FALSE & df_rle2[k, lengths] <= max_gap_num]){
          ind_start = df_rle2[k - 1, cumlength]
          ind_end = df_rle2[k, cumlength] + 1
          if(ind_end > nrow(df_met)) next
          
          df_met[ind_start:ind_end, wnd := approx(datetime, wnd, xout = datetime)$y]
        }
      }
    } # End of gap-filling
    
    # Remove first and last day if they are less than 24 hours
    if(hour(df_exp[1L, datetime]) > 0){
      df_exp = df_exp[datetime >= ceiling_date(df_exp[1L, datetime], unit = "days")]
    }
    if(hour(df_exp[.N, datetime]) <= 23){
      df_exp = df_exp[datetime < floor_date(df_exp[.N, datetime], unit = "days")]
    }
    
    # Compute k600
    df_exp = merge(df_exp, df_met, by = "datetime", all.x = T)
    df_exp[, k600 := k.crusius(df_exp)$k600]
    
    # # Alternatives:
    # df_exp[, k600 := k.cole(df_exp)$k600]
    # df_exp[, k600 := 0.2112 * 1.0241^(Temp - 20)] # Davidson et al 2018 & Liboriussen et al 2010
    # df_exp[is.na(Temp), k600 := 0.2112]
    
    # Compute do.sat
    df_exp[, do_sat := o2.at.sat(df_exp[, .(datetime,
                                            wtr = TW)],
                                 altitude = elev_sites[names(elev_sites) == j])$do.sat]
    
    # Add "mixing depth" (i.e. mesocosm depth)
    df_exp[, zmix := depth_meso]
    
    # Compute metabolism
    if(remove_days_with_na){
      df_exp = df_exp[complete.cases(df_exp)]
    }
    
    df_metab = df_exp[, .(metab = tryCatch({
                                      function_to_fit(do.obs = O2, do.sat = do_sat, k.gas = k600,
                                        z.mix = zmix, irr = PAR, wtr = TW,
                                        datetime = datetime)$metab
                            },
                            error = function(e) as.numeric(c(NA, NA, NA))),
                                                  Mesocosm = unique(Mesocosm)),
                      by = .(floor_date(datetime, unit = "days"),
                             Sensor)]
    setnames(df_metab, old = "floor_date", new = "datetime")
    df_metab[, metab_var := rep(c("GPP", "R", "NEP"), times = .N / 3)]
    df_metab = dcast(df_metab, datetime + Sensor + Mesocosm ~ metab_var, value.var = "metab")
    
    ### Quality control
    if(apply_gpp_r_bounds){
      df_metab[GPP < min_gpp, GPP := NA]
      df_metab[R > max_r, R := NA]
    }
    
    if(!is.null(sensors_to_remove[[j]][[as.character(i)]])){
      sens_to_remove = sensors_to_remove[[j]][[as.character(i)]]
      df_metab[Sensor %in% sens_to_remove, `:=`(GPP = NA,
                                                R = NA,
                                                NEP = NA)]
    }
    
    df_metab[is.na(GPP) | is.na(R), NEP := NA]
    
    ### End of quality control 
    
    # Write output
    fwrite(df_metab, paste0(folder_out, "Metabolism_output_", i, "_", j, suffix, ".csv"))
    
    for(k in c("GPP", "R", "NEP")){
      p1 = ggplot(df_metab) +
        geom_line(aes(datetime, .data[[k]], colour = Mesocosm, group = Sensor)) +
        geom_point(aes(datetime, .data[[k]], colour = Mesocosm, group = Sensor)) +
        labs(title = paste(k, j, i, sep = " - "), x = "date",
             y = paste(k, " (mgO2/l/d)")) +
        theme_classic()
      
      ggsave(paste0(folder_plots, "/Plot_", k, "_", j, "_", i, suffix, ".", fig_type), plot = p1,
             width = fig_width, height = fig_height)
    }
    
    progress = progress + 1
    setTxtProgressBar(progressBar, progress)
  }
}
