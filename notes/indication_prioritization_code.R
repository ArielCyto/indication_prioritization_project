# load env
library(ggpubr)
library(forcats)
library(dplyr)
library(boot, lib.loc = "/usr/lib/R/library")
devtools::load_all("~/capsule/code/cytoreason.ccm.pipeline/") 
source("~/capsule/code/indication_prioritization_project/R/indication_prioritization_func.R") # temporal location for the project functions

# load input data (still need to QC)
# this object has to be on the BQ as a frozen source
gene_exp_data <- get_outputs_dist("wf-3a3109e1f0")$output.rds

# define effector genes by the user
effector_genes <- c('BRCA1','BRCA2', 'HOXB13')

# run results - few examples
rank_groups_by_genes(effector_genes,10)
rank_groups_by_genes(effector_genes,10, 'CI_high')
rank_groups_by_genes('CD4',100)
rank_groups_by_genes('CD4',100, 'median')
rank_groups_by_genes(c('CD4', 'Cd8'),100) # should be failing but suggest alternative names
rank_groups_by_genes(c('CD4', 'CD8A', 'CD8B'),100)
rank_groups_by_genes('BRCA2')
rank_groups_by_genes('not_a_real_gene') # should be failing
rank_groups_by_genes(c('BRCA1','not_a_real_gene','BRCA2'))  # should be failing
