# load env
library(ggpubr)
library(forcats)
library(dplyr)
library(boot, lib.loc = "/usr/lib/R/library")
devtools::load_all("~/capsule/code/cytoreason.ccm.pipeline/") 
source("~/capsule/code/indication_prioritization_project/R/indication_prioritization_func-entrez_id.R") # temporal location for the project functions
unixtools::set.tempdir("~/capsule/scratch/")
# load input data - takes time for the first time
expression_matrix <- read.csv(get_task_outputs("wf-2c26fb32e8",0)["tcga_sample_gene_data.csv"])
annotation_table <-  read.csv(get_task_outputs("wf-2c26fb32e8",0)["tcga_annotation_data.csv"])
colnames(expression_matrix) <- lapply(colnames(expression_matrix), gsub, pattern='X', replacement='')

head(expression_matrix[1:5,1:5])
rownames(expression_matrix) <- annotation_table$sample_id
head(expression_matrix[1:5,1:5])

# fix the numbers for LUAD and LUSC
table(annotation_table$condition)
table(annotation_table$n_sample_condition)
annotation_table$n_sample_condition[annotation_table$condition == "lung cancer - LUAD"] <- "527"
annotation_table$n_sample_condition[annotation_table$condition == "lung cancer - LUSC"] <- "501"
table(annotation_table$n_sample_condition)
# fix the percent sample condition
annotation_table$new_percent_sample_condition <- annotation_table$percent_sample_condition
annotation_table$new_percent_sample_condition <- as.numeric(annotation_table$n_samples)/as.numeric(annotation_table$n_sample_condition)*100
annotation_table$percent_sample_condition <- annotation_table$new_percent_sample_condition
annotation_table$new_percent_sample_condition <- NULL

# add new metadata
# (1) add new metadata column for all of the indications
# (2) for BRCA and CRC change the group label (concat with the new column)

subtypes_table <- read.csv('~/capsule/code/indication_prioritization_project/indication_prioritization_subtypes_table_ver1.csv')
# take only our IO sample_ids
subtypes_table <- subtypes_table[subtypes_table$sample_id %in% annotation_table$sample_id, ]

# organize the table:
temp_annotation_table <- cbind(annotation_table, subtypes_table)
temp_annotation_table$CRC_sub[temp_annotation_table$CRC_sub == 'nan' | temp_annotation_table$CRC_sub == ''] <- "unclassified"
temp_annotation_table$BRCA_sub[temp_annotation_table$BRCA_sub == 'nan' | temp_annotation_table$BRCA_sub == ''] <- "unclassified"
temp_annotation_table$X <- NULL

# check for BRCA and CRC
table(temp_annotation_table$condition, temp_annotation_table$CRC_sub) # we have also MSI and MSS data for bladder and breast
table(temp_annotation_table$condition, temp_annotation_table$BRCA_sub) # we have also data for the other indications
# according to Lital's request - stay only with Basal, Her2, LumA and LumB
temp_annotation_table$BRCA_sub[!temp_annotation_table$BRCA_sub %in% c('Basal', 'Her2', 'LumA', 'LumB')] <- 'unclassified'
table(temp_annotation_table$condition, temp_annotation_table$BRCA_sub)

# remove the information for the relevant indications
# temp_annotation_table$CRC_sub[temp_annotation_table$condition != "colorectal cancer"] <- "unclassified"
# temp_annotation_table$BRCA_sub[temp_annotation_table$condition != "breast cancer"] <- "unclassified"

# give the new names
temp_annotation_table$group_subytpe <- NA
temp_annotation_table$group_subytpe[!temp_annotation_table$condition %in% c('colorectal cancer', 'breast cancer')] <- 'unclassified'
temp_annotation_table$group_subytpe[temp_annotation_table$condition == 'colorectal cancer'] <- temp_annotation_table$CRC_sub[temp_annotation_table$condition == 'colorectal cancer']
temp_annotation_table$group_subytpe[temp_annotation_table$condition == 'breast cancer'] <- temp_annotation_table$BRCA_sub[temp_annotation_table$condition == 'breast cancer']

# generate the new group column
temp_annotation_table$old_group <- temp_annotation_table$group
temp_annotation_table$group <- paste0(temp_annotation_table$old_group,"_",temp_annotation_table$group_subytpe)

# generate new n_samples column
n_samples <- data.frame(table(temp_annotation_table$group))
for (group in temp_annotation_table$group){
  temp_annotation_table$n_samples[temp_annotation_table$group == group] <- n_samples$Freq[n_samples$Var1 == group]
}

# generate new percent_sample_condition column
temp_annotation_table$percent_sample_condition <- as.numeric(temp_annotation_table$n_samples)/as.numeric(temp_annotation_table$n_sample_condition)*100

# remove unnecessary columns
temp_annotation_table$old_group <- NULL
temp_annotation_table$CRC_sub <- NULL
temp_annotation_table$BRCA_sub <- NULL
# temp_annotation_table$group_subytpe <- NULL
head(temp_annotation_table)
temp_annotation_table <- subset(temp_annotation_table, select = -11)

# generate new term_md column
temp_annotation_table$term_id <- paste0(temp_annotation_table$dataset_id,"_",temp_annotation_table$indication,"_",temp_annotation_table$group)


View(data.frame(table(temp_annotation_table$group)))
length(unique(temp_annotation_table$group)) # 60

# remove groups with less than 5 samples
samples_to_remove <- temp_annotation_table$sample_id[temp_annotation_table$n_samples < 5]
temp_annotation_table <- temp_annotation_table[!temp_annotation_table$sample_id %in% samples_to_remove,]
expression_matrix <- expression_matrix[!rownames(expression_matrix) %in% samples_to_remove,]
# we need to stay with 3768 samples instead of 3791 samples


# Upload to Cyto-cc
# Upload the final object to Cyto-cc
tags <- list(list(name="owner", value="Ariel Simon"),
             list(name="project", value="indication prioritation"),
             list(name="data", value="TCGA internal R tables for indication prioritization version 2 - subtypes for BRCA and CRC, filtered for 5 samples"),
             list(name="data_requests", value="https://cytoreason.atlassian.net/wiki/spaces/AR/pages/3705667638/Indication+Prioritization+for+IO+-+R+SDK+tool+-+Engineering+HLD"))

save_as_cyto_cc = function(expression_csv, annotation_csv){
  library(cytoreason.cc.client)
  library(cytoreason.io)
  save_data(expression_csv,"./output/tcga_sample_gene_data.csv", row.names = FALSE)
  save_data(annotation_csv,"./output/tcga_annotation_data.csv", row.names = FALSE)
  
}

# check before send:
head(expression_matrix[1:5,1:5])
temp_annotation_table[1:10,]


wf <- run_function_dist(save_as_cyto_cc,
                        expression_csv=expression_matrix, annotation_csv=temp_annotation_table, tags = tags,
                        image = 'eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:develop_latest') 

# Cyto-CC workflow: wf-cf7fa0cb37
