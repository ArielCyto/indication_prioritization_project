# load env
library(dplyr)
library(tidyr)
devtools::load_all("~/capsule/code/cytoreason.ccm.pipeline/")
# load input data (still need to QC)
# these objects have to be on the BQ as a frozen source
expression_matrix <- get_outputs_dist("wf-219a0c80f1")$output.rds
annotation_table <- get_outputs_dist("wf-c546d82645")$output.rds
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

# Upload the final object to Cyto-cc
tags <- list(list(name="owner", value="Ariel Simon"),
             list(name="project", value="indication prioritation"),
             list(name="data", value="TCGA long format"))

save_as_cyto_cc = function(result_table){
  library(cytoreason.cc.client)
  library(cytoreason.io)
  save_data(result_table,"./output/io_tcga_database.csv")
}

wf <- run_function_dist(save_as_cyto_cc,
                        result_table=result_table, tags = tags,
                        image = 'eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:develop_latest') 

# Cyto-CC workflow: wf-566e268555

cyto_cc_result_table <- read.csv((get_task_outputs("wf-566e268555", "0"))["io_tcga_database.csv"])
  