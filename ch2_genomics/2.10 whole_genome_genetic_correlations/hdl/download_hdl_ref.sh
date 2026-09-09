#!bin/bash
#$ -cwd
#$ -m bea
#$ -l h_rt=1:0:0
#$ -pe smp 1
#$ -l h_vmem=5G
#$ -j y
#$ -N download_hdl_ref
#$ -o download_hdl_ref.log

# Download HDL reference panel

mkdir ~/data/external/genomics/hdl

cd ~/data/external/genomics/hdl

wget -c -t 1 \
https://www.dropbox.com/s/6js1dzy4tkc3gac/UKB_imputed_SVD_eigen99_extraction.tar.gz?dl=0 \
--no-check-certificate -O ./UKB_imputed_SVD_eigen99_extraction.tar.gz

tar -xzvf UKB_imputed_SVD_eigen99_extraction.tar.gz

rm UKB_imputed_SVD_eigen99_extraction.tar.gz
