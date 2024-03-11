library(checkmate)
library(forcats)
library(dplyr)
library(ggpubr)
library(boot, lib.loc = "/usr/lib/R/library")
devtools::load_all("~/capsule/code/cytoreason.ccm.pipeline/") 
# library(cytoreason.analytics)
library(cytoreason.gx)
library(Matrix, lib.loc = "/usr/lib/R/library")
library(checkmate)
source("~/capsule/code/indication_prioritization_project/ssGsea_test/ssGsea_function.R")

# load input data - takes time for the first time
expression_matrix <- read.csv(get_task_outputs("wf-2c26fb32e8",0)["tcga_sample_gene_data.csv"])
annotation_table <-  read.csv(get_task_outputs("wf-2c26fb32e8",0)["tcga_annotation_data.csv"])
colnames(expression_matrix) <- lapply(colnames(expression_matrix), gsub, pattern='X', replacement='')

# define effector genes by the user
effector_genes <- c('59', '72', '8038', '8728', '1264', '1282', '1503', '10272', '3315', '3486', '221749', '8482', '9644', '6876', '7045', '7145', '7168')

# for ssGsea
start.time <- Sys.time()
effector_expression_matrix <- t(expression_matrix) 
effector_genes_for_ssGsea <- list( "effector_genes" = effector_genes)

# fit
res_Ariel <- service_ssgsea(effector_expression_matrix, effector_genes_for_ssGsea)


# combine with the algorithm

# work with smaller object of gene_exp_data with all the necessary columns
data_table <- annotation_table
data_table$effector_genes_ssGsea <- t(res_Ariel)


## BOOTSTRAPING - calculate the CI per type_tme (based on the effector cells average vector)
set.seed(42)

gene_stat <- function(d, i) median(d[i, "effector_genes_ssGsea"]) # median calculation function 
get_statistics <- function(data_table, repeats = 5) {
  # per type_tme in data_table, calculate the "stat" (median) function X repeats and keep the statistical results
  bo <- by(data_table, data_table$group, function(i) boot(i, gene_stat, R = repeats))
  # using the t0 value and other values,
  # calculate confidence interval per type_tme in bo, using boot.ci function
  bo_ci <- t(sapply(bo, function(x) {
    bo_ci <- boot.ci(x, type = "basic")
    # bo_ci keeps 3 fields per type_tme: t0, low_ci and high_ci
    c(t0 = bo_ci$t0, CI_low = bo_ci$basic[4], CI_high = bo_ci$basic[5])
  }))
  bo_ci <- as.data.frame(bo_ci)
  bo_ci$group <- rownames(bo_ci)
  # combine the results with the input data table
  data_table <- left_join(data_table, bo_ci, by = "group")
  return(data_table)
}

data_table <- get_statistics(data_table,repeats =9999)

# organize the results - sum the 4053 samples to 30 rows by the type_tme (change "t0" name to "median")
output_object <-
  data_table %>%
  dplyr::group_by(condition) %>%
  dplyr::group_by(condition, sub_population) %>%
  dplyr::summarise(
    group = unique(group),
    median =unique(t0),
    CI_low = unique(CI_low),
    CI_high = unique(CI_high),
    n_samples = unique(n_samples),
    percent_sample_condition = unique(percent_sample_condition)
    
  )
output_object <- relocate(output_object, group) # put group as the first column

# visualize results
# print results plot

output_object$group = with(output_object, reorder(group, CI_low))
print(output_object %>%
        slice_head(n = 10) %>%
        ggplot(aes(x = CI_low, y = group)) + ggtitle(paste0("Ranking by ", paste(effector_genes, collapse = ', '), " - ssGsea ranking"))+
        geom_point(alpha = 0.8, size = 3, aes(colour = factor(sub_population))) +
        geom_errorbar(aes(xmin = CI_low, xmax = CI_high),
                      width = .2, size = .8,
                      position = "identity", colour = "darkgray"))

# print results numeric (order by "CI_low" as default)
output_object_for_client <- data.frame(output_object[order(output_object$CI_low, decreasing = TRUE),])
output_object_for_client$rank <- 1:(dim(output_object_for_client)[1])
output_object_for_client <- relocate(output_object_for_client, rank)


head(output_object_for_client,nrow(output_object_for_client))
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken
