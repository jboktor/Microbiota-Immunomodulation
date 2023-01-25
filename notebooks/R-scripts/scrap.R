
# First R script
library(reticulate)
# renv::use_python(python= '/home/jboktor/miniconda3/envs/pdmbsR' )

reticulate::use_condaenv(condaenv = 'pdmbsR', required = TRUE)
reticulate::conda_install("pdmbsR", "gget", channel = "bioconda")
