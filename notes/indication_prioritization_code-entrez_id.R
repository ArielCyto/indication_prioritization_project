# load env
library(ggpubr)
library(forcats)
library(dplyr)
library(boot, lib.loc = "/usr/lib/R/library")
devtools::load_all("~/capsule/code/cytoreason.ccm.pipeline/") 
source("~/capsule/code/indication_prioritization_project/R/indication_prioritization_func-entrez_id.R") # temporal location for the project functions

# load input data - takes time for the first time
expression_matrix <- read.csv(get_task_outputs("wf-2c26fb32e8",0)["tcga_sample_gene_data.csv"])
annotation_table <-  read.csv(get_task_outputs("wf-2c26fb32e8",0)["tcga_annotation_data.csv"])
colnames(expression_matrix) <- lapply(colnames(expression_matrix), gsub, pattern='X', replacement='')
# define effector genes by the user
effector_genes <- c("1234", "3122", "10563", "6352", "6772", "4818", "4283", "3903", "3902", "10663", "3824",
                    "3001", "5551", "100049587", "26191", "942", "6503", "55423", "971", "84868", "9050", "114836",
                    "8832", "146722", "915", "3458", "6373", "914", "1522", "3002", "3561", "3627", "11006", "5133",
                    "6355", "4261", "6351", "10261", "5788", "283420", "8530", "3620", "3683", "999", "9051",
                    "3003", "3133", "916", "117289", "3604")

# run results - few examples
stromal_genes <- c("2191",'59', '7040', '5154','1277', '1278', '1281', '4313', '6387', '4318',
           '2246','6347','6935', '910' ,"1301", "3624" ,"7058" ,"1290")
rank_groups_by_genes(stromal_genes,9999)
