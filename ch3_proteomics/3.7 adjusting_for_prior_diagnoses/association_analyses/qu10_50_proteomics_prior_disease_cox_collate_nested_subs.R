library(tidyverse)

disease <- as.character(commandArgs(trailingOnly = TRUE)[3])

protein_result_dir <- paste0('~/ch3_proteomics/3.7 adjusting_for_prior_diagnoses/association_analyses/cox_results/', disease)

# this should be done by counting the number of files - because Firth error runs still export an empty csv, and the last protein may not be the last to finish
protein_result_files <- list.files(protein_result_dir, full.names = TRUE, recursive = FALSE)
if (length(protein_result_files) == 2920){
  all_cox_results <- data.frame()
  for (file in list.files(protein_result_dir, full.names = TRUE)){
    all_cox_results <- rbind(all_cox_results, read.csv(file))
  }
  write_csv(all_cox_results, paste0('~/ch3_proteomics/3.7 adjusting_for_prior_diagnoses/association_analyses/cox_results/', disease, '_all_proteins.csv'))
  
  print(paste0(Sys.time(), ': Analyses complete. Results and fitted models exported.'))
  flush.console()
  Sys.sleep(1)
  
}

