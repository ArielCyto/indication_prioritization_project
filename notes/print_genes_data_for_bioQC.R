# load env
library(ggpubr)
library(forcats)
library(dplyr)
library(boot, lib.loc = "/usr/lib/R/library")
devtools::load_all("~/capsule/code/cytoreason.ccm.pipeline/") 
source("~/capsule/code/indication_prioritization_project/R/indication_prioritization_func.R") # temporal location for the project functions

# load input data (still need to QC)
# these objects have to be on the BQ as a frozen source
expression_matrix <- get_outputs_dist("wf-219a0c80f1")$output.rds
annotation_table <- get_outputs_dist("wf-c546d82645")$output.rds

for (disease_name in unique(annotation_table$indication)){
  disease_samples <- (annotation_table[annotation_table$indication == disease_name,])$sample_id
  disease_popular_genes <- expression_matrix[rownames(expression_matrix) %in% disease_samples,]
  disease_top_genes <- sort(colMeans(disease_popular_genes),decreasing = TRUE)
  print(paste0("The top 20 genes for ",disease_name, " are:"))
  print(data.frame(head(disease_top_genes, 20)))
}


# check random gene list
random_genes <- sample(x=colnames(expression_matrix), size=19)
rank_groups_by_genes(random_genes)
