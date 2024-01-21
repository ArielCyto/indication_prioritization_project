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

# Checking median / average with Lital
stromal_genes <- c('FAP','ACTA2','TGFB1','PDGFA','COL1A1','COL1A2','COL3A1','MMP2','CXCL12','MMP9','FGF1','CCL2','ZEB1','CD1B','COL11A1','INHBA','THBS2','COL5A2')
rank_groups_by_genes(stromal_genes)
inflamed_markers <- c('PDCD1','IFNG','STAT1','CCL5','CXCL9','CXCL10','HLA-A','TNFAIP1','IDO1','PDCD1LG2','STAT3','FOXP3','CD4','CD3D','CD8A','CD274','CTLA4','LAG3','ACOX2')
rank_groups_by_genes(inflamed_markers)
lital_test1 <- c('SIRPA','CD47', 'CD80', 'CD86', 'CD163')
rank_groups_by_genes(lital_test1)
lital_test2 <- c('CD28', 'CTLA4', 'CD27', 'PDCD1', 'CBLB', 'EGF', 'CD47', 'TGFB1')
rank_groups_by_genes(lital_test2)

