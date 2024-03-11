# load env
library(ggpubr)
library(forcats)
library(dplyr)
library(boot, lib.loc = "/usr/lib/R/library")
remotes::install_github("s-u/unixtools")
unixtools::set.tempdir("/scratch")
devtools::load_all("~/capsule/code/cytoreason.ccm.pipeline/") 
source("~/capsule/code/indication_prioritization_project/R/indication_prioritization_func-entrez_id.R") # temporal location for the project functions


# load input data - takes time for the first time
# expression_matrix <- read.csv(get_task_outputs("wf-2c26fb32e8",0)["tcga_sample_gene_data.csv"])
annotation_table <-  read.csv(get_task_outputs("wf-2c26fb32e8",0)["tcga_annotation_data.csv"])
TCGA_updated <- readRDS(get_workflow_outputs('wf-01d457dac0'))
expression_matrix <- TCGA_updated@assayData$counts # TPM
head(expression_matrix[1:5,1:5])

# colnames(expression_matrix) <- lapply(colnames(expression_matrix), gsub, pattern='X', replacement='')

# find the entrez id of GZMA and PRF1
# library(org.Hs.eg.db)
# hs <- org.Hs.eg.db
# my.symbols <- c("GZMA", "PRF1")
# select(hs, 
#        keys = my.symbols,
#        columns = c("ENTREZID", "SYMBOL"),
#        keytype = "SYMBOL")
# # SYMBOL ENTREZID
# # 1   GZMA     3001
# # 2   PRF1     5551

# to avoid negative or 0 values, we will add the minimal number in the table + 0.00001 to have only positive numbers
expression_matrix <- expression_matrix + abs(min(expression_matrix)) + 0.00001

# calc the geometric mean per sample
GZMA_PRF1 <- c('3001', '5551')
GZMA_PRF1_vals <- expression_matrix[,GZMA_PRF1]
GZMA_PRF1_vals$geometric_mean <- GZMA_PRF1_vals$`3001` # just to initialize the new column
for(i in 1:nrow(GZMA_PRF1_vals)) {
  GZMA_PRF1_vals$geometric_mean[i] <- exp(mean(log(as.numeric(GZMA_PRF1_vals[i,1:2]))))
}

# normalize all the matrix by the geometric mean of the sample 
expression_matrix <- expression_matrix/GZMA_PRF1_vals$geometric_mean

# bioQC for SO usecase - genes of mature natural killer cell
IO_signatures <- readRDS(get_workflow_outputs("wf-d2ad9d9a50"))
sign <- IO_signatures@.Data
names(sign)<- IO_signatures@names
IO_signatures <- sign

genes_to_check <- sign$`mature natural killer cell`
genes_to_check <- genes_to_check[genes_to_check %in% colnames(expression_matrix)]
dfi <- rank_groups_by_genes(genes_to_check)

genes_to_check <- sign$`macrophage`
genes_to_check <- genes_to_check[genes_to_check %in% colnames(expression_matrix)]
rank_groups_by_genes(genes_to_check)

