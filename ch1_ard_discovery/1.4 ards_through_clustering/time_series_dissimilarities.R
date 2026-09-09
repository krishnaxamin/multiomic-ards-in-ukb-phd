# install.packages('tidyverse')
# install.packages('TSclust')

library(tidyverse)
library(TSclust)

' Script calculating dissimilarity between all time series in a given directory 
using 7 dissimiarlity measures from the R package TSclust. '

# dissimilarities calculated using the TSclust package
# diss_measure taken from 'distances' vector, e.g. eucl, not Euclidean
# which_time_series = 'original', 'smoothed', 'normalised_smoothed'
# time_series_path must have '/' at the end
dissimilarity_calc <- function(diss_measure, which_time_series = 'normalised_smoothed', time_series_path = '~/data/internal/time_series/genomics_qcv1_time_series/raw'){
  
  time_series_file_names <- as.list(list.files(path = time_series_path, pattern = '*_time_series.csv'))
  num_time_series <- length(time_series_file_names)
  num_time_series_pairs <- choose(num_time_series, 2)
  
  stats <- rep(0, num_time_series_pairs)
  stat_counter <- 1
  
  for(i in seq_along(l=time_series_file_names)) {
    if (i != num_time_series){
      for (file_name in time_series_file_names[(i+1):num_time_series]){
        time_series_pair <- c(time_series_file_names[[i]], file_name[[1]])
        time_series_1 <- read.csv(paste('time_series/', time_series_pair[1], sep = ''))
        time_series_2 <- read.csv(paste('time_series/', time_series_pair[2], sep = ''))
        time_series_variants <- list('original' = list(x = time_series_1$disease_onset_rate,
                                                       y = time_series_2$disease_onset_rate),
                                     'smoothed' = list(x = time_series_1$smoothed_profile,
                                                       y = time_series_2$smoothed_profile),
                                     'normalised_smoothed' = list(x = time_series_1$normalised_smoothed,
                                                                  y = time_series_2$normalised_smoothed))
        time_series_x = time_series_variants[[which_time_series]]$x
        time_series_y = time_series_variants[[which_time_series]]$y
        
        if (diss_measure == 'cdm'){
          stat <- diss.CDM(time_series_x, 
                           time_series_y)
        }
        else if (diss_measure == 'cid'){
          stat <- diss.CID(time_series_x, 
                           time_series_y)
        }
        else if (diss_measure == 'cort'){
          stat <- diss.CORT(time_series_x, 
                            time_series_y,
                            k = 2,
                            deltamethod = 'Euclid')
        }
        else if (diss_measure == 'dtw'){
          stat <- diss.DTWARP(time_series_x, 
                              time_series_y)
        }
        else if (diss_measure == 'eucl'){
          stat <- diss.EUCL(time_series_x, 
                            time_series_y)
        }
        else if (diss_measure == 'frechet'){
          stat <- diss.FRECHET(time_series_x, 
                               time_series_y)
        }
        else if (diss_measure == 'pearson'){
          stat <- diss.COR(time_series_x, 
                           time_series_y)
        }
        else {
          print('Dissimiliarty measure not one of: cdm, cid, cort, dtw, eucl, frechet, lpc, pearson.')
        }
        stats[stat_counter] <- stat
        stat_counter = stat_counter + 1
      }
    }
  }
  
  # save diss stats
  saveRDS(stats, paste('~/ch1_ard_discovery/1.4 ards_through_clustering/distances/normalised_smoothed_', diss_measure, 'cort_distances.rds', sep = ''))
  
  print(paste(diss_measure, 'complete.'))
  
}

distances = c('cdm', 'cid', 'cort', 'dtw', 'eucl', 'frechet', 'pearson')

for (distance in distances){
  dissimilarity_calc(distance)
}
