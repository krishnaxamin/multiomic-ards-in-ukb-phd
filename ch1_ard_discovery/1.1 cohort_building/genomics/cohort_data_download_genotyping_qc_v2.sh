#!bin/bash
#cohort_data_download for Genotyping QC v2
# use Table Exporter app to download a table of desired data fields from a cohort
# the desired data fields are listed in a unix-format text file (field_names_file_txt) - one field per line - array-genotyping_phenotype_fields_to_download.txt
dx run table-exporter \
-idataset_or_cohort_or_dashboard='{"$dnanexus_link": "record-..."}' \
-ientity=participant \
-ifield_names_file_txt='{"$dnanexus_link": "file-..."}' \
-ioutput=array-genotyping_phenotype_file \
-y --instance-type mem3_ssd1_v2_x8 --priority low --name array-genotyping_phenotype_extraction