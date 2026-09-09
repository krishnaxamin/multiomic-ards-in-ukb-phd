#!bin/bash
#$ -cwd
#% -m bea
#$ -l h_rt=1:0:0
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l centos
#$ -j y
#$ -N magma_step1_annotate
#$ -o ~/gwas/magma/mamga_step1_annotate.log

# annotate SNPs to genes using MAGMA for MAGMA gene-level analysis

# navigate to directory where the MAGMA executable is
# cd ~/gwas/magma
./magma \
--annotate window=5 \
--snp-loc ~/data/external/genomics/magma/g1000_eur.bim \
--gene-loc ~/data/external/genomics/NCBI37.3.gene.loc \
--out ~/data/external/magma/all_genotyped_imputed
