# Scripts to create a sensor data.table in a consistent format, to unify the
# different formats from the different years and sites
# This script uses the 1-minute (2-minute for Skogaryd) files

library(data.table)
library(stringr)

Sys.setenv(TZ = "UTC")

source("o2sat_to_conc.R")

compile_sensor_ts = function(folder_data, site, year){
  # Find filename
  abbr = fcase(site == "Erken", "ERK",
               site == "Bolmen", "BOL",
               site == "Skogaryd", "SRC")
  filename = list.files(folder_data, pattern = paste0("SITES_PHY_", abbr))
  filename = filename[grepl(year, filename)]
  filename = filename[grepl("min\\.csv", filename)]
  
  lines_file = readLines(file.path(folder_data, filename))
  skiprows = which(lines_file == "####")
  df = fread(file.path(folder_data, filename), skip = skiprows)
  df[, TIMESTAMP := as.POSIXct(TIMESTAMP, format = "%Y-%m-%d %H:%M")]
  
  # In 2023 files, have to calculate O2 concentration from saturation
  if(sum(grepl("O2_", names(df))) == 0L){
    o2sat_cols = names(df)[grepl("O2SAT_", names(df))]
    for(k in o2sat_cols){
      temp_col = gsub("O2SAT_", "TW_", k)
      new_col = gsub("O2SAT_", "O2_", k)
      df[, (new_col) := o2sat_to_conc(get(k), get(temp_col))]
    }
  }
  
  # Keep relevant columns
  cols_to_keep = grepl("TIMESTAMP|PAR|O2|TW", names(df))
  df = df[, ..cols_to_keep]
  
  # Split into treatments and different variables
  df_long = melt(df, id.vars = "TIMESTAMP", variable.factor = F,
                 variable.name = "Treatment")
  df_long[, Variable := sapply(str_split(Treatment, "_"), "[[", 1)]
  df_long[, Sensor := fifelse(grepl("Lake", Treatment), "lake", "mesocosm")]
  df_long[Sensor == "mesocosm", Sensor := str_match(Treatment, "_(\\d+)_")[, 2]]
  df_long[Sensor == "lake", Mesocosm := "Lake"]
  df_long[Sensor != "lake", Mesocosm := str_match(Treatment, "_(\\w)$")[, 2]]
  
  df_long[, Mesocosm := fcase(Mesocosm == "Lake", "Lake",
                              Mesocosm == "C", "Control",
                              Mesocosm == "D", "Daily",
                              Mesocosm == "I", "Intermediate",
                              Mesocosm == "E", "Extreme")]
  
  df_exp = dcast(df_long, TIMESTAMP + Sensor + Mesocosm ~ Variable,
                 value.var = "value", fun.aggregate = mean)
  setnames(df_exp, old = "TIMESTAMP", "datetime")
  
  return(df_exp)
}
