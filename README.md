# predator-prey_RCdates
Data and code to the manuscript "A note on predator prey dynamics in radiocarbon datasets"

INSTRUCTIONS

to run the analysis, 

1) make sure that the following libraries are installed: 

rcarbon
philentropy
infotheo
modelbased


2) Load the data object "20231109_data.Rdata" into the global environment. The object contains vectors of radiocarbon dates and errors by trophic level and location, as well as the "desco_trophic" and "fairbanks_trophic" databases. 

3) Run the function kld_dates() in 20231101_code.R. The header of the file contains ab explanation of the function paramenters. 

4) To reproduce the Fairbanks results, run

kld_dates(fairbanks_prey$RC_date, fairbanks_prey$RC_error, fairbanks_predator$RC_date, fairbanks_predator$RC_error, sample = 100, dataset_name = "Fairbanks", y_scaling_parameter = 4, smoothing = 0.4)

5) To reproduce the Judean Desert results, run 

kld_dates(desco_prey$RC_date, desco_prey$RC_error, desco_predator$RC_date, desco_predator$RC_error, sample = 100, dataset_name = "Judean Desert", y_scaling_parameter = 4, smoothing = 0.4)
