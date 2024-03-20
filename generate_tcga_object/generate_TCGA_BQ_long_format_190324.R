# load env
library(dplyr)
library(tidyr)
devtools::load_all("~/capsule/code/cytoreason.ccm.pipeline/")
install.packages('unixtools', repos = 'http://www.rforge.net/')
library(unixtools)
unixtools::set.tempdir("/scratch")

# these objects have to be on the BQ as a frozen source
expression_matrix <- read.csv(get_task_outputs("wf-cf7fa0cb37",0)["tcga_sample_gene_data.csv"])
annotation_table <-  read.csv(get_task_outputs("wf-cf7fa0cb37",0)["tcga_annotation_data.csv"])
colnames(expression_matrix) <- lapply(colnames(expression_matrix), gsub, pattern='X', replacement='')


table(annotation_table$group)

# change column names according to https://cytoreason.atlassian.net/wiki/spaces/AR/pages/3705667638/Indication+Prioritization+for+IO+-+R+SDK+tool+-+Engineering+HLD
# Reshape the data to long format
rownames(expression_matrix) <- annotation_table$sample_id
expression_matrix$sample_id <- rownames(expression_matrix)
head(expression_matrix[1:5,1:5])

head(annotation_table) # change the bad column typo :) !!!!
colnames(annotation_table)[11] <- "group_subtype"

# Multiple variables stored in column names
long_data <- expression_matrix %>% pivot_longer(
  cols = -sample_id,
  names_to = "feature_id",
  values_to = "values"
)

# change the name values to value
head(long_data)
colnames(long_data)[3] <- "value"
head(long_data)

# Add columns: domain and dataset_id
long_data$domain <- "io"
long_data$dataset <- "TCGA"
head(long_data)

# add groups
# Merge with metadata
groups_data <- data.frame(sample_id = annotation_table$sample_id, group = annotation_table$group)
merged_data <- left_join(long_data,groups_data , by = "sample_id")
head(merged_data)

# check NA
isTRUE(is.na(merged_data)) # need to return false





# annotation table:
annotation_table_for_cyto_cc <- subset(annotation_table, select = -sample_id)
rownames(annotation_table_for_cyto_cc) <- NULL
annotation_table_for_cyto_cc <- distinct(annotation_table_for_cyto_cc)

head(annotation_table_for_cyto_cc)

# Upload the final object to Cyto-cc
tags <- list(list(name="owner", value="Ariel Simon"),
             list(name="project", value="indication prioritation"),
             list(name="data", value="TCGA final tables for indication prioritation version 2 - long format, fix typo"),
             list(name="data_requests", value="https://cytoreason.atlassian.net/wiki/spaces/AR/pages/3705667638/Indication+Prioritization+for+IO+-+R+SDK+tool+-+Engineering+HLD"))

save_as_cyto_cc = function(expression_csv, annotation_csv){
  library(cytoreason.cc.client)
  library(cytoreason.io)
  save_data(expression_csv,"./output/tcga_sample_gene_data.csv", row.names = FALSE)
  save_data(annotation_csv,"./output/tcga_annotation_data.csv", row.names = FALSE)
  
}

wf <- run_function_dist(save_as_cyto_cc,
                        expression_csv=merged_data, annotation_csv=annotation_table_for_cyto_cc, tags = tags,
                        image = 'eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:develop_latest') 
# Cyto-CC workflow: wf-83a0c1c25e - with the typo
# Cyto-CC workflow: wf-bc9435edd8