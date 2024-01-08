# load env
library(ggpubr)
library(forcats)
library(dplyr)
library(boot, lib.loc = "/usr/lib/R/library")
devtools::load_all("~/capsule/code/cytoreason.ccm.pipeline/") 

# load input data
gene_exp_data <- get_outputs_dist("wf-3a3109e1f0")$output.rds

# define effector genes by the user
effector_genes <- c('BRCA1','BRCA2', 'HOXB13')



# check if the genes are available:
submit = FALSE
for (gene in effector_genes){
  gene_status <- gene %in% colnames(gene_exp_data)
  if (gene_status){
    submit = TRUE
  }else{
    print(paste0("gene ",gene," is an invalid input, please remove it and try again" ))
    submit = FALSE
    # check if you can suggest alternative names for the wrong gene name
    optional_alternative_name <- colnames(gene_exp_data)[(grepl(toupper(gene),colnames(gene_exp_data)))]
    if (length(optional_alternative_name) > 0){
      print("optional alternatives: ")
      print(optional_alternative_name)
    }
    
    break
  }
}

# run analysis
if (submit == TRUE){
  
  # per sample (row), calculate the selected effector_genes average expression
  if (length(effector_genes) > 1) # multiple gene list
  {gene_exp_data$effector_genes_average <- rowMeans(gene_exp_data[,effector_genes])
  }else{ # one gene only
    gene_exp_data$effector_genes_average <- gene_exp_data[,effector_genes]
  }
  
  
  # work with smaller object of gene_exp_data with all the necessary columns
  data_table <- data.frame(type_tme = gene_exp_data$type_tme, cluster_annotation = gene_exp_data$cluster_annotation,
                           tissue_io = gene_exp_data$tissue_io,
                           effector_genes_average = gene_exp_data$effector_genes_average,
                           n_per_indication_TME = gene_exp_data$n_per_indication_TME,
                           percent_tme = gene_exp_data$percent_tme)
  
  
  ## BOOTSTRAPING - calculate the CI per type_tme (based on the effector cells average vector)
  set.seed(42)
  
  gene_stat <- function(d, i) median(d[i, "effector_genes_average"]) # median calculation function 
  get_statistics <- function(data_table, repeats = 5) {
    # per type_tme in data_table, calculate the "stat" (median) function X repeats and keep the statistical results
    bo <- by(data_table, data_table$type_tme, function(i) boot(i, gene_stat, R = repeats))
    # using the t0 value and other values,
    # calculate confidence interval per type_tme in bo, using boot.ci function
    bo_ci <- t(sapply(bo, function(x) {
      bo_ci <- boot.ci(x, type = "basic")
      # bo_ci keeps 3 fields per type_tme: t0, low_ci and high_ci
      c(t0 = bo_ci$t0, CI_low = bo_ci$basic[4], CI_high = bo_ci$basic[5])
    }))
    bo_ci <- as.data.frame(bo_ci)
    bo_ci$type_tme <- rownames(bo_ci)
    # combine the results with the input data table
    data_table <- left_join(data_table, bo_ci, by = "type_tme")
    return(data_table)
  }
  
  data_table <- get_statistics(data_table,repeats =20)
  
# organize the results - sum the 4053 samples to 30 rows by the type_tme (change "t0" name to "score")
  output_object <-
    data_table %>%
    dplyr::group_by(tissue_io) %>%
    dplyr::group_by(tissue_io, cluster_annotation) %>%
    dplyr::summarise(
      type_tme = unique(type_tme),
      score =unique(t0),
      CI_low = unique(CI_low),
      CI_high = unique(CI_high),
      n_per_indication_TME = unique(n_per_indication_TME),
      percent_tme = unique(percent_tme)

    )
  output_object <- relocate(output_object, type_tme) # put type_tme as the first column
  
  # visualize results
  # print results plot
  
  output_object$type_tme = with(output_object, reorder(type_tme, score))
  print(output_object %>%
    slice_head(n = 10) %>%
    ggplot(aes(x = score, y = type_tme)) + ggtitle(paste0("Ranking by ", paste(effector_genes, collapse = ', ')))+
    geom_point(alpha = 0.8, size = 3, aes(colour = factor(cluster_annotation))) +
    geom_errorbar(aes(xmin = CI_low, xmax = CI_high),
                  width = .2, size = .8,
                  position = "identity", colour = "darkgray"))
  
  # print results numeric
  output_object_for_client <- data.frame(output_object[order(output_object$score, decreasing = TRUE),])
  head(output_object_for_client,30)
  
}


