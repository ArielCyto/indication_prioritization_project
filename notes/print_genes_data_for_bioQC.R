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
annotation_table <- get_outputs_dist("wf-1e6d943964")$output.rds # with the LUNG splitted !

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
Immune_check_point <- c('CD274','PDCD1LG2','CTLA4','PDCD1','LAG3','HAVCR2','TIGIT')
rank_groups_by_genes(Immune_check_point)
MKI67_cell_cycle <- c('MKI67', 'CCNE1', 'BUB1', 'BUB1B', 'CCNB2', 'CDC25C', 'CDK2', 'MCM4', 'MCM6',
                      'MCM2')
rank_groups_by_genes(MKI67_cell_cycle)
TGFB_receptor_ligand <- c('TGFB1', 'TGFBR2')
rank_groups_by_genes(TGFB_receptor_ligand)
FTBRS <- c('ACTA2','COL4A1', 'TAGLN', 'SH3PXD2A')
rank_groups_by_genes(FTBRS)

# When entrez_id are asked:
gene_to_entrez <- read.csv("capsule/code/indication_prioritization_project/gene-to-entrez.csv")
gene_to_entrez <-  subset(gene_to_entrez, select = -X)
result_table$entrez_id = result_table$gene_id # generate the new column

convert_entrez_to_genes <- function(entrez_list){
  gene_list <- vector(mode='list', length=length(entrez_list))
  counter = 1
  for (entrez in entrez_list){
    index = gene_to_entrez[gene_to_entrez$entrez_id == entrez,]
    if (dim(index)[1] != 0){
      gene_list[counter] <- as.character(index$gene_name)
      counter = counter+1
    }
    else{
      print(paste0(gene," is problematic"))
    }
  }
  
  return(unlist(gene_list))
}


IFNG_list <- c("1234", "3122", "10563", "6352", "6772", "4818", "4283", "3903", "3902", "10663", "3824", "3001", "5551", "100049587", "26191", "942", "6503", "55423", "971", "84868", "9050", "114836", "8832", "146722", "915", "3458", "6373", "914", "1522", "3002", "3561", "3627", "11006", "5133", "6355", "4261", "6351", "10261", "5788", "283420", "8530", "3620", "3683", "999", "9051", "3003", "3133", "916", "117289", "3604")
rank_groups_by_genes(convert_entrez_to_genes(IFNG_list))

TGFB_list <- c('59', '72', '8038', '8728', '1264', '1282', '1503', '10272', '3315', '3486', '221749', '8482', '9644', '6876', '7045', '7145', '7168')
rank_groups_by_genes(convert_entrez_to_genes(TGFB_list))
