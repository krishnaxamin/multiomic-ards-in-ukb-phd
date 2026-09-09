#!bin/bash
# Download metabolomics data for all those with the data. data covers covariates, metabolite data, QC flags
# use Table Exporter app to download a table of desired data fields from a cohort
# the desired data fields are listed in a unix-format text file (field_names_file_txt) - one field per line - nmr_metabolomics_fields_to_download_no_empty_qc_fields.txt
dx run table-exporter \
-idataset_or_cohort_or_dashboard='{"$dnanexus_link": "record-..."}' \
-ientity=participant \
-ifield_names_file_txt='{"$dnanexus_link": "file-..."}' \
-ioutput=metabolomics_data_covars_qcflags_off_dnanexus \
-y --instance-type mem3_ssd1_v2_x8 --priority low --name nmr_metabolomics_data_extraction