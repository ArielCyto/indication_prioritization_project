# load env
library(dplyr)
library(tidyr)
devtools::load_all("~/capsule/code/cytoreason.ccm.pipeline/")
install.packages('unixtools', repos = 'http://www.rforge.net/')
library(unixtools)
unixtools::set.tempdir("/scratch")

# these objects have to be on the BQ as a frozen source
expression_matrix <- get_outputs_dist("wf-219a0c80f1")$output.rds
expression_matrix <- subset(expression_matrix, select = -c(effector_genes_average))
annotation_table <- get_outputs_dist("wf-1e6d943964")$output.rds

# Reshape the data to long format
expression_matrix$group <-  annotation_table$group

# Add entrez
gene_to_entrez <- read.csv("capsule/code/indication_prioritization_project/gene-to-entrez.csv")
gene_to_entrez <-  subset(gene_to_entrez, select = -X)
gene_to_entrez$gene_name <- as.character(gene_to_entrez$gene_name)
gene_to_entrez[nrow(gene_to_entrez)+1,] = c("no gene","nothing")

expression_matrix[nrow(expression_matrix) + 1,] = gene_to_entrez$entrez_id # takes 1 min
rownames(expression_matrix[nrow(expression_matrix),]) <- "entrez_id"
expression_matrix$sample_id <- rownames(expression_matrix)


# long_data <- pivot_longer(expression_matrix, cols = -sample_id, names_to = c("gene_id", 'entrez_id', 'group'), values_to = "values")

# Multiple variables stored in column names
long_data <- expression_matrix %>% pivot_longer(
  cols = -sample_id,
  names_to = "gene_id",
  values_to = "values"
)

# rm(expression_matrix); gc() # cus the memory is out of bond
# Merge with metadata
colnames(gene_to_entrez)[1] <- 'gene_id'
merged_data <- left_join(long_data, gene_to_entrez, by = "gene_id")
groups_data <- data.frame(sample_id = annotation_table$sample_id, group = annotation_table$group)
merged_data <- left_join(merged_data,groups_data , by = "sample_id")
head(merged_data)

# remove ".1" for dup genes c('ABCF2.1', 'ABCF2', 'ECE2', 'ECE2.1', 'HSPA14', 'HSPA14.1', 'POLR2J3', 'POLR2J3.1',
# 'TMSB15B', 'TMSB15B.1')
merged_data$gene_id[merged_data$gene_id == "ABCF2.1"] <- "ABCF2"
merged_data$gene_id[merged_data$gene_id == "ECE2.1"] <- "ECE2"
merged_data$gene_id[merged_data$gene_id == "HSPA14.1"] <- "HSPA14"
merged_data$gene_id[merged_data$gene_id == "POLR2J3.1"] <- "POLR2J3"
merged_data$gene_id[merged_data$gene_id == "TMSB15B.1"] <- "TMSB15B"


# Upload the final object to Cyto-cc
tags <- list(list(name="owner", value="Ariel Simon"),
             list(name="project", value="indication prioritation"),
             list(name="data", value="TCGA final table1"))

save_as_cyto_cc = function(result_table){
  library(cytoreason.cc.client)
  library(cytoreason.io)
  save_data(result_table,"./output/tcga_sample_gene_data.csv", row.names = FALSE)
}

wf <- run_function_dist(save_as_cyto_cc,
                        result_table=merged_data, tags = tags,
                        image = 'eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:develop_latest') 

# Cyto-CC workflow: wf-21b6e14f82
# https://cyto-cc.cytoreason.com/workflow/wf-21b6e14f82


# annotation table:
annotation_table_for_cyto_cc <- subset(annotation_table, select = -sample_id)
rownames(annotation_table_for_cyto_cc) <- NULL
annotation_table_for_cyto_cc <- distinct(annotation_table_for_cyto_cc)

# Upload the final object to Cyto-cc
tags <- list(list(name="owner", value="Ariel Simon"),
             list(name="project", value="indication prioritation"),
             list(name="data", value="TCGA final table2"))

save_as_cyto_cc = function(result_table){
  library(cytoreason.cc.client)
  library(cytoreason.io)
  save_data(result_table,"./output/tcga_annotation_data.csv", row.names = FALSE)
}

wf <- run_function_dist(save_as_cyto_cc,
                        result_table=annotation_table_for_cyto_cc, tags = tags,
                        image = 'eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:develop_latest') 