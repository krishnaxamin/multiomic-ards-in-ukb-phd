#!bin/bash
# Download proteomics data for all those with the data and that are sex concordant
# use Table Exporter app to download a table of desired data fields from a cohort
# the desired data fields are listed in a unix-format text file (field_names_file_txt) - one field per line - proteomics_olinks_fields_to_download.csv & proteomics_participant_fields_to_download

# First, download the proteomics data from the olink_instance_0 entity
dx run table-exporter \
-idataset_or_cohort_or_dashboard='{"$dnanexus_link": "record-..."}' \
-ientity=olink_instance_0 \
-ifield_names_file_txt='{"$dnanexus_link": "file-..."}' \
-ioutput=proteomics_npx_data \
-y --instance-type mem3_ssd1_v2_x8 --priority low --name proteomics_npx_data_extraction

# Then, download the non-proteomics data that contributes to covariates from the participant entity
dx run table-exporter \
-idataset_or_cohort_or_dashboard='{"$dnanexus_link": "record-..."}' \
-ientity=participant \
-ifield_names_file_txt='{"$dnanexus_link": "file-..."}' \
-ioutput=proteomics_non_npx_data \
-y --instance-type mem3_ssd1_v2_x8 --priority low --name proteomics_non_npx_data_extraction