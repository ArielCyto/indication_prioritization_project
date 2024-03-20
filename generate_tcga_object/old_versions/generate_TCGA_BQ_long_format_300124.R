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
table(annotation_table$group) # check the LUAD, LUSC

# Add entrez
gene_to_entrez <- read.csv("capsule/code/indication_prioritization_project/gene-to-entrez.csv")
# gene_to_entrez <-  subset(gene_to_entrez, select = -X)
colnames(expression_matrix) <- gene_to_entrez$entrez_id
# check 
expression_matrix[1:5,1:5]


# add new TME names
annotation_table$TME <- as.character(annotation_table$TME)
annotation_table$TME[annotation_table$TME == "stromal epithelial"] <- "stromal_epithelial"
annotation_table$TME[annotation_table$TME == "dom lymph myeloid"] <- "inflamed"
annotation_table$TME[annotation_table$TME == "dom lymph myeloid mild stromal"] <- "inflamed_mild_stromal"
annotation_table$TME[annotation_table$TME == "stromal dom myeloid lymph"] <- "inflamed_stromal"
table(annotation_table$TME) # check

# change the corresponding group names
annotation_table$group <- paste0(annotation_table$indication,"_",annotation_table$TME)
table(annotation_table$group)

# Reshape the data to long format
expression_matrix$sample_id <- rownames(expression_matrix)


# Multiple variables stored in column names
long_data <- expression_matrix %>% pivot_longer(
  cols = -sample_id,
  names_to = "entrez_id",
  values_to = "values"
)

# change the name values to value
head(long_data)
colnames(long_data)[3] <- "value"
head(long_data)

# add groups


# rm(expression_matrix); gc() # cus the memory is out of bond
# Merge with metadata
groups_data <- data.frame(sample_id = annotation_table$sample_id, group = annotation_table$group)
merged_data <- left_join(long_data,groups_data , by = "sample_id")
head(merged_data)

# check NA
isTRUE(is.na(merged_data)) # need to false

# annotation table:
annotation_table_for_cyto_cc <- subset(annotation_table, select = -sample_id)
rownames(annotation_table_for_cyto_cc) <- NULL
annotation_table_for_cyto_cc <- distinct(annotation_table_for_cyto_cc)
colnames(annotation_table_for_cyto_cc)[6] <- "percent_TME_in_indication"
  
# Upload the final object to Cyto-cc
tags <- list(list(name="owner", value="Ariel Simon"),
             list(name="project", value="indication prioritation"),
             list(name="data", value="TCGA final tables for indication prioritation"))

save_as_cyto_cc = function(expression_csv, annotation_csv){
  library(cytoreason.cc.client)
  library(cytoreason.io)
  save_data(expression_csv,"./output/tcga_sample_gene_data.csv", row.names = FALSE)
  save_data(annotation_csv,"./output/tcga_annotation_data.csv", row.names = FALSE)
  
}

wf <- run_function_dist(save_as_cyto_cc,
                        expression_csv=merged_data, annotation_csv=annotation_table_for_cyto_cc, tags = tags,
                        image = 'eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:develop_latest') 
# https://cyto-cc.cytoreason.com/workflow/wf-3dfbbd2ea9