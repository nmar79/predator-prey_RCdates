#Calculates the KL(Predator || Prey), when Predator and Prey are radiocarbon date summed probability distributions (SPDs). The function also provides bootstrap support in relation to randomly produced SPDs from the same time range of Prey. 
#prey_dates: first sample of uncalibrated 14C dates. This is the reference distribution.
#predator_dates: second sample of uncalibrated radiocarbon dates
#prey_errors, predator_errors: associated 14C measurement errors. These errors also serve to estimate errors for random datasets used to produce the bootstrapped SPD sample.
#n_samples: number of random samples of N = length(predator_dates) to construct random SPDs simulating Predator. 
#dataset_name: optional title for the resulting plots. Defaults to "".
#smoothing: loess smoothing intensity, from 0 to 1. Calls modelbased::smoothing. Defaults to 0. 
#y_scaling_parameter: scaling factor for the y-axis. Must be larger than 0, defaults to 1. 


kld_dates <- function(prey_dates, prey_errors, predator_dates, predator_errors, dataset_name = "", sample_n = 100, smoothing = 0, y_scaling_parameter = 1){

require(rcarbon)
require(philentropy)
require(infotheo)
require(modelbased)
  
if (y_scaling_parameter <= 0){stop("y_scaling_parameter must be larger than 0")}
  
if (sample_n >1000){stop("Too many iterations, aborting...")}
  
if (sample_n > 100){
  
  invisible(readline(prompt="Many iterations, it can be slow. Press [Enter] to continue"))
 
} 
 
set.seed(42)
par(mfrow = c(2,1))

#determine range

earliest_date <- max(prey_dates, predator_dates)
latest_date <- min(prey_dates, predator_dates)

  #calibrates and extracts the values of the prey SPD

  print("Calibrating prey dates and producing SPD...")  
  calib_prey <- calibrate(x = prey_dates, errors = prey_errors, verbose = FALSE)
  prey_spd <- spd(calib_prey, timeRange = c(earliest_date, latest_date), spdnormalised = TRUE, verbose = FALSE)
  prey_spd_values <- prey_spd$grid$PrDens
  
  print("Done.")
  
  #calibrates and extracts the values of the predator SPD
  
  print("Calibrating predator dates and producing SPD...")
  
  calib_predator <- calibrate(x = predator_dates, errors = predator_errors, verbose = FALSE)
  predator_spd <- spd(calib_predator, timeRange = c(earliest_date, latest_date), spdnormalised = TRUE, verbose = FALSE)
  predator_spd_values <- predator_spd$grid$PrDens
  
  print("Done.")
  
  #stash max density for plot
  y_high <- max(predator_spd_values, prey_spd_values)
  
  #get Kullback-Leibler Divergence statistic between the prey and predator SPDs 
  
  print("calculating Kullback-Leibler Divergence between prey and predator SPDs...")
  
  origin_kl_mx <- rbind(prey_spd_values, predator_spd_values)
  origin_kl <- suppressMessages(KL(origin_kl_mx))
  
  print("Done.")
  


#The "produce_error" function takes a random date, finds the closest archaeological wolf date, and returns the error of that wolf date. 

print("loading 'produce_error' helper function to calculate likely errors to be associated with random dates...")
  
produce_error <- function(rand_date, predator_dates){
  dist_rand_real <- 1:length(predator_dates)
  for (i in 1:length(predator_dates)){dist_rand_real[i] <- sqrt((rand_date - predator_dates[i])^2)}
  closest_error <- min(dist_rand_real)
  temp_counter <- which(dist_rand_real == closest_error)
  rand_error <- predator_errors[temp_counter]
  return(rand_error)
}

print("Done.")

#associate random date with semi-random error for a random sample

print("Sampling...")


kl <- 1:sample_n

progress_bar <- txtProgressBar(min = 0, max = sample_n, style = 3, width = 50, char = "=")   

for (i in 1:sample_n){

sample_dates <- 1:length(predator_dates)
sample_errors <- 1:length(predator_dates)

for (j in 1:length(predator_dates)){

rand_date <- sample(latest_date:earliest_date, 1, replace = TRUE)

while(rand_date < 0){rand_date <- sample(latest_date:earliest_date, 1, replace = TRUE)}

rand_error <- produce_error(rand_date, predator_dates)

sample_dates[j] <- rand_date
suppressWarnings(sample_errors[j] <- rand_error)
}

###diagnostics: plot(sample_dates, sample_errors) -- clear!

#calibrate the N=10 sample and calculate spd

calib_sample <- calibrate(x = sample_dates, errors = sample_errors, verbose = FALSE)
sample_spd <- spd(calib_sample, timeRange = c(earliest_date,latest_date), spdnormalised = TRUE, verbose = FALSE)

if(i == 1){
plot(y=modelbased::smoothing(sample_spd$grid$PrDens, strength = smoothing), x = prey_spd$grid$calBP, ylim = c(0, y_high/y_scaling_parameter),type = "l", xlab = "years cal BP", ylab = "probability density", col = "lightgrey", lwd = 0.5, main = dataset_name)} 
else 
  {lines(y=modelbased::smoothing(sample_spd$grid$PrDens, strength = smoothing), x = prey_spd$grid$calBP, type = "l", col = "lightgrey", lwd = 0.5)}


#COMPARISON
prey_pd <- prey_spd_values
sample_pd <- sample_spd$grid$PrDens

KL_mx <- rbind(prey_pd, sample_pd)
kl[i] <- suppressMessages(KL(KL_mx))

setTxtProgressBar(progress_bar, i)

}

close(progress_bar)

print("Plot prey and predator SPDs...")

lines(modelbased::smoothing(prey_spd_values, strength = smoothing), x = prey_spd$grid$calBP, type = "l", lwd = 2, col = "darkgreen")
lines(y = modelbased::smoothing(predator_spd_values,strength = smoothing), x = prey_spd$grid$calBP, col = "darkred", lwd = 2)

print("Done.")

print("_____________RESULTS_________________")

print(paste0("Summary statistics for Kullback-Leibler Divergence"," ","(original distance = ", round(origin_kl,4), ")"))

print(summary(kl))

x <- which(kl > origin_kl)
p <- length(x)/length(kl)


print(paste0("The Kullback-Leibler Divergence between the original datasets is smaller than between the original curves in ", p ," of the replicates"))
  
hist(kl, main = "", xlab = "Kullback-Leibler Divergence")
abline(v=origin_kl, col = "red")

}

