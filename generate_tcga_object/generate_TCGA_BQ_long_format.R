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
expression_matrix$sample_id <- rownames(expression_matrix)
long_data <- pivot_longer(expression_matrix, cols = -sample_id, names_to = "gene_id", values_to = "values")
# rm(expression_matrix); gc() # cus the memory is out of bond
# Merge with metadata
colnames(annotation_table)[colnames(annotation_table) == '%_TME_in_indication'] <- 'percent_TME_in_indication'
merged_data <- left_join(long_data, annotation_table, by = "sample_id")

# Group by gene_id and the annotations
result_table <- merged_data %>%
  group_by(gene_id, TME, indication, group, n_per_indication, n_samples, percent_TME_in_indication) %>%
  summarize(
    sample_ids = toString(sample_id),
    values = toString(values)
  )


########## Try adding entrez_id
gene_to_entrez <- read.csv("capsule/code/indication_prioritization_project/gene-to-entrez.csv")
gene_to_entrez <-  subset(gene_to_entrez, select = -X)
result_table$entrez_id = result_table$gene_id # generate the new column


for (gene in unique(result_table$gene_id)){
  index = gene_to_entrez[gene_to_entrez$gene_name == gene,]
  if (dim(index)[1] != 0){
    result_table$entrez_id[result_table$gene_id == gene] <- index$entrez_id
  }
  else{
    print(paste0(gene," is problematic"))
  }
}

# remove ".1" for dup genes c('ABCF2.1', 'ABCF2', 'ECE2', 'ECE2.1', 'HSPA14', 'HSPA14.1', 'POLR2J3', 'POLR2J3.1',
# 'TMSB15B', 'TMSB15B.1')
# result_table$gene_id[result_table$gene_id == "ABCF2.1"] <- "ABCF2"
# result_table$gene_id[result_table$gene_id == "ECE2.1"] <- "ECE2"
# result_table$gene_id[result_table$gene_id == "HSPA14.1"] <- "HSPA14"
# result_table$gene_id[result_table$gene_id == "POLR2J3.1"] <- "POLR2J3"
# result_table$gene_id[result_table$gene_id == "TMSB15B.1"] <- "TMSB15B"


# Upload the final object to Cyto-cc
tags <- list(list(name="owner", value="Ariel Simon"),
             list(name="project", value="indication prioritation"),
             list(name="data", value="TCGA long format - LUAD / LUSC splitted, entrez_id added"))

save_as_cyto_cc = function(result_table){
  library(cytoreason.cc.client)
  library(cytoreason.io)
  save_data(result_table,"./output/io_tcga_database.csv", row.names = FALSE)
}

wf <- run_function_dist(save_as_cyto_cc,
                        result_table=result_table, tags = tags,
                        image = 'eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:develop_latest') 

# Cyto-CC workflow:wf-1d9e387afe
# https://cyto-cc.cytoreason.com/workflow/wf-1d9e387afe