# Coxmos install

# install Bioconductor for this version of R, e.g. version 3.9 for R4.4
install.packages('BiocManager')
BiocManager::install(version = "3.9")

# install mixOmics - do not update packages
BiocManager::install('mixOmics')

# install survcomp - do not update packages
BiocManager::install('survcomp')

# install Coxmos
install.packages('Coxmos')