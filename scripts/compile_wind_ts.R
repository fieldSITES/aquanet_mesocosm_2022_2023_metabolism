# Scripts to create a wind data.table in a consistent format, to unify the
# different formats from the different years and sites

library(data.table)
library(readxl)

Sys.setenv(TZ = "UTC")

compile_wind_ts = function(folder_data, site, year){
  if(!(year %in% c(2022, 2023))){
    stop("'year' should be either '2022' or '2023'!")
  }
  
  if(site == "Erken"){
    if(year == 2022){
      df = fread(file.path(folder_data, "Meteo_Erken_2022.dat"))
      
      df = df[, .(datetime = TMSTAMP,
                  wnd = MeanWS)]
    }else if(year == 2023){
      df = fread(file.path(folder_data, "Erken_AirTemp_BaroPressure_Humidity_2023.csv"))
      
      df = df[, .(datetime = TIMESTAMP,
                  wnd = MeanWS)]
    }
  }else if(site == "Bolmen"){
    df_ws = read_excel(file.path(folder_data, "smhi-opendata_3_4_63510_20240527_045035.xlsx"),
                       skip = 10L)
    setDT(df_ws)
    df_ws[, datetime := as.POSIXct(paste0(Datum, " ", hour(`Tid (UTC)`), ":",
                                          minute(`Tid (UTC)`), ":", second(`Tid (UTC)`)))]
    df = df_ws[, .(datetime,
                   wnd = as.numeric(Vindhastighet))]
  }else if(site == "Skogaryd"){
    if(year == 2022){
      stop("No experiment was done in 2022 in Skogaryd!")
    }
    df = read_excel(file.path(folder_data, "Skogaryd_Raft_wsense_aquacosm.xlsx"))
    setDT(df)
    df = df[, .(datetime = TIMESTAMP,
                wnd = WS_avg)]
  }else{
    stop("'site' argument should be one of 'Erken', 'Bolmen', or 'Skogaryd'!")
  }
  
  # Subset time
  if(year == 2022){
    df = df[datetime >= as.POSIXct("2022-07-01") & datetime <= as.POSIXct("2022-08-31")]
  }else if(year == 2023){
    df = df[datetime >= as.POSIXct("2023-04-01") & datetime <= as.POSIXct("2023-06-30")]
  }
  
  return(df)
}
