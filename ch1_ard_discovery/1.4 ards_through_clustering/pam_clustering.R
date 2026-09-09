library(tidyverse)
library(cluster) # for PAM and clusGap

' Code to perform clustering, given an optimal k and a dissimilarity measure. 
Code also to plot the clusters and further information about the clusters. '

# CORT ------

cort_stats <- readRDS('~/ch1_ard_discovery/1.4 ards_through_clustering/distances/normalised_smoothed_cort_distances.rds')

# PAM clustering
clustering <- pam(cort_stats, k = 4, diss = TRUE)

# note to which cluster each disease belongs
disease_in_which_cluster_df <- data.frame(disease = disease_names,
                                          cluster = clustering$clustering)
colnames(disease_in_which_cluster_df) <- c('disease_field', 'cluster')

# include ICD-10 information
disease_cluster_icd10_df <- merge(disease_in_which_cluster_df, 
                                  disease_fields_to_icd10,
                                  by = 'disease_field')

write.csv(disease_cluster_icd10_df, 
          '~/ch1_ard_discovery/1.4 ards_through_clustering/clusters/disease_in_which_cluster_cort_pam-k4_norm_smooth_no-harmonics.csv',
          row.names=FALSE)

## Plotting ------

time_series_file_names <- as.list(list.files(path = '~/data/internal/time_series/genomics_qcv1_time_series/raw/', pattern = '*_time_series.csv'))
num_time_series <- length(time_series_file_names)
disease_in_which_cluster_df <- read.csv('~/ch1_ard_discovery/1.4 ards_through_clustering/clusters/disease_in_which_cluster_cort_pam-k4_norm_smooth_no-harmonics.csv')
disease_field_chapters <- read.csv('~/data/external/phenotype_coding/disease_fields_icd10_chapters.csv')

# set up character vector as a lookup table to give the cluster plots labels describing how many diseases are in each cluster
# number of diseases in each cluster
num_diseases_per_cluster <- table(disease_in_which_cluster_df$cluster)
cluster_plot_labels <- rep(0, length(num_diseases_per_cluster))
for (i in 1:length(num_diseases_per_cluster)) {
  i_ch <- as.character(i)
  cluster_plot_labels[i] <- (i_ch = paste('Cluster ', i_ch, ': ', num_diseases_per_cluster[i], ' diseases', sep = ''))
}

df_to_plot <- data.frame(age = numeric(),
                         disease_onset_rate = numeric(),
                         disease = character(),
                         cluster = integer(),
                         icd10_chapter = integer())

# num_chapter_in_cluster <- rep(0, num_time_series)

for(i in seq_along(l=time_series_file_names)) {
  time_series <- read.csv(paste('~/data/internal/time_series/genomics_qcv1_time_series/raw/', time_series_file_names[[i]], sep=''))
  disease_field <- str_split(time_series_file_names[[i]], '_')[[1]][1]
  cluster <- disease_in_which_cluster_df[disease_in_which_cluster_df$disease == disease_field, "cluster"]
  cluster_label <- cluster_plot_labels[cluster]
  icd10_chapter <- disease_field_chapters[disease_field_chapters$disease_field == disease_field, "icd10_chapter"]
  
  time_series$disease <- disease_field
  time_series$cluster <- cluster
  time_series$icd10_chapter <- icd10_chapter
  time_series$cluster_label <- cluster_label
  
  df_to_plot <- rbind(df_to_plot, time_series)
  # num_chapter_in_cluster[i] <- paste('chapter ', as.character(icd10_chapter), ' in cluster ', as.character(cluster), sep='')
}

df_to_plot$cluster <- as.factor(df_to_plot$cluster)
df_to_plot$cluster_label <- as.factor(df_to_plot$cluster_label)
df_to_plot$icd10_chapter <- as.factor(df_to_plot$icd10_chapter)

ggplot(data = df_to_plot, aes(x = age, y = normalised_smoothed, group = disease, color = icd10_chapter)) +
  geom_line() + 
  facet_wrap(vars(cluster_label)) +
  labs(x = 'Age',
       y = 'Disease onset rate',
       title = 'Diseases clustered by normalised smoothed age-of-onset profile (harmonics removed)',
       subtitle = 'Dissimilarity measure: CORT. Clustering algo: PAM, k=4.') + 
  scale_colour_discrete(name = 'ICD10 chapter')
ggsave(paste('~/ch1_ard_discovery/1.4 ards_through_clustering/plots/normalised-smooth-profiles/cort_pam/cort_pam-k4_clusters.pdf', sep = ''),
       device = 'pdf', width = 297, height = 210, units = 'mm')
ggsave(paste('~/ch1_ard_discovery/1.4 ards_through_clustering/plots/normalised-smooth-profiles/cort_pam/cort_pam-k4_clusters.png', sep = ''),
       device = 'png', width = 297, height = 210, units = 'mm')
