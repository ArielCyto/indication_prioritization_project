# load env
library(ggpubr)
library(forcats)
library(dplyr)
library(boot, lib.loc = "/usr/lib/R/library")
devtools::load_all("~/capsule/code/cytoreason.ccm.pipeline/") 


# load the TCGA genes object 
annotation_table <- read.csv("capsule/code/indication_prioritization_project/annotation_table_wf_2441e859da.csv")
expression_matrix <- readRDS("capsule/code/indication_prioritization_project/expression_matrix_wf_e692424459.rds")

# remove uncommon sample_ids
samples_to_drop <- annotation_table$sample_id[!annotation_table$sample_id %in% rownames(expression_matrix)]
annotation_table[annotation_table$sample_id %in% samples_to_drop,]
annotation_table <- annotation_table[!annotation_table$sample_id %in% samples_to_drop,]
# check dim
dim(annotation_table)[1] == dim(expression_matrix)[1]

# change to SYMBOL - I already changed it 
head(colnames(expression_matrix))

# generate unexisted columns
annotation_table$group <- paste0(annotation_table$indication,"_",annotation_table$TME)
table(annotation_table$group)
length(unique(annotation_table$group))
# [1] 26

# Add columns: n_per_indication, n_samples, %_TME_in_indication
annotation_table <- data.frame(annotation_table)
annotation_table <- annotation_table %>%
  # generate n_per_indication and n_samples columns
  dplyr::add_count(indication, name = "n_per_indication") %>%
  dplyr::add_count(indication, TME, name = "n_samples") %>%
  # generate %_TME_in_indication by the two variables above
  dplyr::mutate('%_TME_in_indication' = (n_samples / n_per_indication) * 100)

annotation_table$X <- NULL
rownames(annotation_table) <- annotation_table$sample_id


# Upload the final objects to cyto-cc
tags <- list(list(name="owner", value="Ariel Simon"),
             list(name="project", value="IO indication prioritization"),
             list(name="dataset", value="TCGA"),
             list(name="data", value="annotation_table")
             )

wf <- run_function_dist(function(obj){return(obj)},
                        obj=annotation_table, tags = tags,
                        image = 'eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:develop_latest') 
# Cyto-CC workflow: wf-c546d82645
  

tags <- list(list(name="owner", value="Ariel Simon"),
             list(name="project", value="IO indication prioritization"),
             list(name="dataset", value="TCGA"),
             list(name="data", value="expression_matrix")
)

wf <- run_function_dist(function(obj){return(obj)},
                        obj=expression_matrix, tags = tags,
                        image = 'eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:develop_latest') 
# Cyto-CC workflow: wf-219a0c80f1


