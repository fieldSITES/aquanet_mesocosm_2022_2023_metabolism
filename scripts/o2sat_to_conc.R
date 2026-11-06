# Function to convert oxygen saturation to concentration, form the Aanderaa
# oxygen sensor manual

o2sat_to_conc = function(o2sat, temp){
  ts = log((298.15 - temp)/ (273.15 + temp))
  
  a0 = 2.00856
  a1 = 3.22400
  a2 = 3.99063
  a3 = 4.80299
  a4 = 0.978188
  a5 = 1.71069
  
  cstar = exp(a0 + a1 * ts + a2 * ts^2 + a3 * ts^3 + a4 * ts^4 + a5 * ts^5)
  
  doconc_uM = cstar * 44.659 * o2sat / 100
  
  return(doconc_uM * 32 / 1000) # conversion to mg/l
}
