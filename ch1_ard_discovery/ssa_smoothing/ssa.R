# install.packages('Rssa')
# install.packages('tidyverse')

library(Rssa)
library(tidyverse)

source('~/ch1_ard_discovery/1.4 ards_through_clustering/ssa_smoothing/ssa_periodicity_detection_func.R')

' Script exploring single spectrum analysis. A function for SSA on one time 
series is defined and carried out on the full set of time series. '


# SSA pipeline function (defined for one time series) ------

ssa_pipeline <- function(time_series_for_func){
  
  # produce ssa object
  ssa_test <- ssa(time_series_for_func$disease_onset_rate, kind = '1d-ssa', svd.method = 'auto')
  
  # get percentage-variance-explained by each eigenvector
  eigenvals <- (ssa_test$sigma)^2
  percentage_variance_explained <- eigenvals / sum(eigenvals)
  
  # choose eigenvectors based on their percentage-variance-explained
  desired_Fn <- which(percentage_variance_explained > 0.01)
  
  # filter chosen eigenvectors by their periodicity
  desired_Fn <- intersect(desired_Fn, which(apply(ssa_test$U, MARGIN = 2, gradient_switch_counter) < 10))
  
  # get how much variance is explained by the chosen eigenvectors
  percentage_variance_explained_by_chosen_Fn <- sum(percentage_variance_explained[desired_Fn])
  
  # number of chosen eigenvectors
  num_Fn <- length(desired_Fn)
  
  reconstruction_pipeline <- reconstruct(ssa_test, group = list(Trend = desired_Fn))
  # plot(reconstruction_pipeline, add.residuals = FALSE, main = 'Pipeline test')
  
  smoothed_profile <- reconstruction_pipeline$Trend 
  smoothed_profile[smoothed_profile < 0] <- 0 # converts negative values to zeroes - mull over benefit of doing this
  
  time_series_for_func$smoothed_profile <- smoothed_profile
  
  return_object <- list(time_series = time_series_for_func, 
                        total_var_explained = percentage_variance_explained_by_chosen_Fn,
                        num_eigenvecs_chosen = num_Fn)
  return_object
}

# SSA Pipeline (for many time series) ------

## Raw time series ------

time_series_file_names <- as.list(list.files(path = '~/data/internal/time_series/genomics_qcv1_time_series/raw/', pattern = '*_time_series.csv'))
num_time_series <- length(time_series_file_names)

disease_names <- rep(0, num_time_series)
all_total_var_explained <- rep(0, num_time_series)
all_num_Fn <- rep(0, num_time_series)

for(i in seq_along(l=time_series_file_names)) {
  time_series <- read.csv(paste('~/data/internal/time_series/genomics_qcv1_time_series/raw/', time_series_file_names[[i]], sep = ''))
  disease_names[i] <- str_split(time_series_file_names[[i]], '_')[[1]][1]
  ssa_result <- ssa_pipeline(time_series)
  time_series_with_smoothed <- ssa_result$time_series 
  
  # normalise so smoothed_profile sums to 1
  time_series_with_smoothed <- mutate(time_series_with_smoothed, normalised_smoothed = smoothed_profile / sum(smoothed_profile))
  
  write.csv(time_series_with_smoothed, file = paste('~/data/internal/time_series/genomics_qcv1_time_series/raw/', time_series_file_names[[i]], sep = ''), row.names = FALSE)
  all_total_var_explained[i] <- ssa_result$total_var_explained
  all_num_Fn[i] <- ssa_result$num_eigenvecs_chosen
}
