# load env
library(ggpubr)
library(forcats)
library(dplyr)
library(boot, lib.loc = "/usr/lib/R/library")
devtools::load_all("~/capsule/code/cytoreason.ccm.pipeline/") 
source("~/capsule/code/indication_prioritization_project/R/indication_prioritization_func-entrez_id.R")
unixtools::set.tempdir("~/capsule/scratch/")

# load input data - the first time takes time
expression_matrix <- read.csv(get_task_outputs("wf-cf7fa0cb37",0)["tcga_sample_gene_data.csv"])
annotation_table <-  read.csv(get_task_outputs("wf-cf7fa0cb37",0)["tcga_annotation_data.csv"])
colnames(expression_matrix) <- lapply(colnames(expression_matrix), gsub, pattern='X', replacement='')
rownames(expression_matrix) <- annotation_table$sample_id

# example 1 - check for stromal genes list
effector_genes <- c("2191",'59', '7040', '5154','1277', '1278', '1281', '4313', '6387', '4318',
                   '2246','6347','6935', '910' ,"1301", "3624" ,"7058" ,"1290")
rank_groups_by_genes(effector_genes,9999) # 9999 is the number of repeats in the bootstrap process.


# example 2 - check for 1 gene and ask for less repeats (only for out internal use to save time)
rank_groups_by_genes('7040',5)
