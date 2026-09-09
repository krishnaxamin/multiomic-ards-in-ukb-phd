library(tidyverse)
library(cluster) # for PAM and clusGap
library(factoextra) # for fviz_gap_stat
library(ggpubr) # for ggarrange

' Script to find the optimal k for the PAM algorithm. '


optimal_pam_k <- function(distance, distance_name){
  
  library('factoextra')
  library('ggpubr')
  
  distance_stats = readRDS(paste('~/ch1_ard_discovery/1.4 ards_through_clustering/distances/normalised_smoothed_', distance, '_distances.rds', sep = ''))
  
  # build, swap are the values of the objective function of the PAM algo after the first two steps. Equiv. to WCC for k-means clustering.
  build_log <- rep(0, 10)
  swap_log <- rep(0, 10)
  
  # look at average silhouette width over a range of k
  sil_avg_width_log <- rep(0, 10)
  
  # collect plots looking at silhouette widths more in detail
  plot_list <- list()
  
  for (i in 1:10){
    # print(i)
    set.seed(42)
    pam_result <- pam(distance_stats, k = i, diss = TRUE)
    build_log[i] <- pam_result$objective['build'][[1]]
    swap_log[i] <- pam_result$objective['swap'][[1]]
    if (is.null(pam_result$silinfo) == FALSE) {
      sil_avg_width_log[i] <- pam_result$silinfo$avg.width[[1]]
      silhouette_plot_title <- paste(distance_name, '. PAM. k = ', i, '. Avg. silhouette width: ', round(pam_result$silinfo$avg.width[[1]], 3), sep = '')
      plot_list[[i]] <- fviz_silhouette(pam_result, print.summary = FALSE, legend = 'none', title = silhouette_plot_title)
    } else {
      sil_avg_width_log[i] <- 0
    }
  }
  
  # plot the objective functions
  pam_diff_ks <- data.frame(k = c(1:10),
                            build = build_log,
                            swap = swap_log) %>% 
    pivot_longer(!k, names_to = 'measure', values_to = 'value')
  pam_diff_ks$measure <- as.factor(pam_diff_ks$measure)
  ggplot(pam_diff_ks, aes(x = factor(k, level = 1:10), y = value, group = measure, colour = measure)) + 
    geom_line() + 
    geom_point() +
    labs(x = 'k',
         y = 'Value',
         title = 'Objective function values over k',
         subtitle = paste('Dissimilarity measure: ', distance_name, '. Clustering algo: PAM.', sep = '')) +
    scale_colour_discrete(name = 'Measure')
  ggsave(paste('~/ch1_ard_discovery/1.4 ards_through_clustering/plots/normalised-smooth-profiles/', distance, '_pam/', distance, '_pam_obj-fun-over-k.pdf', sep = ''),
         device = 'pdf', width = 297, height = 210, units = 'mm')
  ggsave(paste('~/ch1_ard_discovery/1.4 ards_through_clustering/plots/normalised-smooth-profiles/', distance, '_pam/', distance, '_pam_obj-fun-over-k.png', sep = ''),
         device = 'png', width = 297, height = 210, units = 'mm')
  
  print(paste(distance_name, 'objective functions done.'))
  
  # plot the avg silhouette widths
  pam_avg_sil_width <- data.frame(k = c(1:10),
                                  avg_sil_width = sil_avg_width_log)
  ggplot(pam_avg_sil_width, aes(x = factor(k, level = 1:10), y = avg_sil_width)) + 
    geom_line(colour = 'red', group = 1) + 
    geom_point(colour = 'red') + 
    labs(x = 'k', 
         y = 'Average silhouette width',
         title = 'Average silhouette width over k',
         subtitle = paste('Dissimilarity measure: ', distance_name, '. Clustering algo: PAM.', sep = ''))
  ggsave(paste('~/ch1_ard_discovery/1.4 ards_through_clustering/plots/normalised-smooth-profiles/', distance, '_pam/', distance, '_pam_avg-sil-width-over-k.pdf', sep = ''),
         device = 'pdf', width = 297, height = 210, units = 'mm')
  ggsave(paste('~/ch1_ard_discovery/1.4 ards_through_clustering/plots/normalised-smooth-profiles/', distance, '_pam/', distance, '_pam_avg-sil-width-over-k.png', sep = ''),
         device = 'png', width = 297, height = 210, units = 'mm')
  
  print(paste(distance_name, 'avg sil widths done.'))
  
  saveRDS(plot_list, paste('~/ch1_ard_discovery/1.4 ards_through_clustering/plots/normalised-smooth-profiles/', distance, '_pam/', distance, '_pam_sil-widths-over-k.rds', sep = ''))
  
  # plot the silhouette widths - can add a title but has to be done somewhat manually (https://community.rstudio.com/t/adding-main-title-using-ggarrange/85629/2)
  ggarrange(plot_list[[2]],
            plot_list[[3]],
            plot_list[[4]],
            plot_list[[5]],
            plot_list[[6]],
            plot_list[[7]],
            plot_list[[8]],
            plot_list[[9]],
            plot_list[[10]],
            ncol = 2, nrow = 5)
  ggsave(paste('~/ch1_ard_discovery/1.4 ards_through_clustering/plots/normalised-smooth-profiles/', distance, '_pam/', distance, '_pam_sil-widths-over-k.pdf', sep = ''),
         device = 'pdf', width = 297, height = 210, units = 'mm')
  ggsave(paste('~/ch1_ard_discovery/1.4 ards_through_clustering/plots/normalised-smooth-profiles/', distance, '_pam/', distance, '_pam_sil-widths-over-k.png', sep = ''),
         device = 'png', width = 297, height = 210, units = 'mm')
  
  print(paste(distance_name, 'sil widths done.'))
  
  # gap stat
  
  stat_matrix <- matrix(0, nrow = num_time_series, ncol = num_time_series)
  
  vector_index <- 0
  for (i in 1:(num_time_series-1)) {
    # print(paste('i =', i, sep=''))
    for (j in (i+1):num_time_series) {
      # print(paste('j =', j, sep=''))
      vector_index <- vector_index + 1
      stat_matrix[i, j] <- distance_stats[vector_index]
      stat_matrix[j, i] <- distance_stats[vector_index]
    }
  }
  
  pam_for_clusgap <- function(x,k) list(cluster = pam(x,k, cluster.only = TRUE))
  
  clusgap_b1000_power2 <- clusGap(stat_matrix, FUN = pam_for_clusgap, K.max = 10, B = 1000, d.power = 2)
  
  saveRDS(clusgap_b1000_power2, paste('~/ch1_ard_discovery/1.4 ards_through_clustering/gap_stats/normalised-smoothed-profiles/', distance, '_pam_b1000_power2.rds', sep = ''))
  
  clusgap_b1000_power2_plot <- fviz_gap_stat(clusgap_b1000_power2, maxSE = list(method = 'Tibs2001SEmax', SE.factor = 1)) +
    ggtitle(paste('Dissimilarity measure: ', distance_name, '. Clustering algo: PAM. B = 1000, d.power = 2', sep = ''))
  clusgap_b1000_power2_plot
  ggsave(paste('~/ch1_ard_discovery/1.4 ards_through_clustering/plots/normalised-smooth-profiles/', distance, '_pam/', distance, '_pam_gap-stats-over-k.pdf', sep = ''),
         device = 'pdf', width = 297, height = 210, units = 'mm')
  ggsave(paste('~/ch1_ard_discovery/1.4 ards_through_clustering/plots/normalised-smooth-profiles/', distance, '_pam/', distance, '_pam_gap-stats-over-k.png', sep = ''),
         device = 'png', width = 297, height = 210, units = 'mm')
  
  # find all k for which gap(k) >= gap(k+1) - s(k+1). K.max = 10
  print(paste('Dissimilarity measure: ', distance_name, '. Clustering algo: PAM. B = 1000, d.power = 2', sep = ''))
  for (i in 1:9) {
    gap_i <- clusgap_b1000_power2_plot$data$gap[i]
    gap_i1 <- clusgap_b1000_power2_plot$data$gap[i+1]
    s_i1 <- clusgap_b1000_power2_plot$data$SE.sim[i+1]
    tibshirani_satisfied <- gap_i >= gap_i1 - s_i1
    if (tibshirani_satisfied == TRUE) {
      print(paste('k = ', i, ': Tibshirani inequality satisfied.', sep = ''))
    }
  }
  
  print(paste(distance_name, 'gap stat done.'))
  
}

distances = c('cdm', 'cid', 'cort', 'dtw', 'eucl', 'frechet', 'pearson')
distances_names = c('CDM', 'CID', 'CORT', 'DTW', 'Euclidean', 'Frechet', 'Pearson')

for (d in 1:length(distances)){
  optimal_pam_k(distances[d], distances_names[d])
}