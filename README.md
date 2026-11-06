# Metabolic rate calculations for SITES Aquanet mesocosm experiment 2022/2023
The R scripts on this repository support the metabolism estimates (GPP/R/NEP) in Langenheder et al. ("Lake ecosystem responses to variability in dissolved organic matter and nutrient loading across time and space")

Please first install the necessary packages via the "Install_packages.R" script in the `scripts/` folder. The main script to run is the `LakeMetabolizer modelling.R` script. It uses the `LakeMetabolizer` R package (Winslow et al., 2016).

The required sensor data files can be downloaded from the SITES data portal: https://doi.org/10.23700/JSX1-XK65. The scripts assume that these are stored in the `SensorData`folder. 

The script generates csv files with the calculated metabolic rates and several plots. Note that the plots in Langenheder et al. are different, but based on the same data.

## References
Winslow, L.A., Zwart, J.A., Batt, R.D., Dugan, H.A., Woolway, R.I., Corman, J.R., & Read, J.S. (2016). LakeMetabolizer: An R package for estimating lake metabolism from free-water oxygen using diverse statistical models. <i>Inland Waters</i>, 6(4), doi:10.1080/IW-6.4.883
