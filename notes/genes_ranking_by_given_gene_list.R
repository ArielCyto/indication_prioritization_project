# load env
library(ggpubr)
library(forcats)
library(dplyr)
library(boot, lib.loc = "/usr/lib/R/library")
devtools::load_all("~/capsule/code/cytoreason.ccm.pipeline/") 

# load input data
expression_matrix <- get_outputs_dist("wf-219a0c80f1")$output.rds
annotation_table <- get_outputs_dist("wf-c546d82645")$output.rds

# define effector genes by the user
effector_genes <- c('ABCA1',	'ABCA2', 'ABCA3', 'NAT2')


# check if the genes are available:
submit = FALSE
for (gene in effector_genes){
  gene_status <- gene %in% colnames(expression_matrix)
  if (gene_status){
    submit = TRUE
  }else{
    print(paste0("gene ",gene," is an invalid input, please remove it and try again" ))
    submit = FALSE
    # check if you can suggest alternative names for the wrong gene name
    optional_alternative_name <- colnames(expression_matrix)[(grepl(toupper(gene),colnames(expression_matrix)))]
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
  {expression_matrix$effector_genes_average <- rowMeans(expression_matrix[,effector_genes])
  }else{ # one gene only
    expression_matrix$effector_genes_average <- expression_matrix[,effector_genes]
  }
  
  
  # work with smaller object of gene_exp_data with all the necessary columns
  data_table <- annotation_table
  data_table$effector_genes_average <- expression_matrix$effector_genes_average
  
  
  ## BOOTSTRAPING - calculate the CI per type_tme (based on the effector cells average vector)
  set.seed(42)
  
  gene_stat <- function(d, i) median(d[i, "effector_genes_average"]) # median calculation function 
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
  
  data_table <- get_statistics(data_table,repeats =20)
  
# organize the results - sum the 4053 samples to 30 rows by the type_tme (change "t0" name to "median")
  output_object <-
    data_table %>%
    dplyr::group_by(indication) %>%
    dplyr::group_by(indication, TME) %>%
    dplyr::summarise(
      group = unique(group),
      median =unique(t0),
      CI_low = unique(CI_low),
      CI_high = unique(CI_high),
      n_samples = unique(n_samples),
      '%_TME_in_indication' = unique(`%_TME_in_indication`)

    )
  output_object <- relocate(output_object, group) # put group as the first column
  
  # visualize results
  # print results plot
  
  output_object$group = with(output_object, reorder(group, CI_low))
  print(output_object %>%
    slice_head(n = 10) %>%
    ggplot(aes(x = CI_low, y = group)) + ggtitle(paste0("Ranking by ", paste(effector_genes, collapse = ', ')))+
    geom_point(alpha = 0.8, size = 3, aes(colour = factor(TME))) +
    geom_errorbar(aes(xmin = CI_low, xmax = CI_high),
                  width = .2, size = .8,
                  position = "identity", colour = "darkgray"))
  
  # print results numeric (order by "CI_low" as default)
  output_object_for_client <- data.frame(output_object[order(output_object$CI_low, decreasing = TRUE),])
  output_object_for_client$rank <- 1:(dim(output_object_for_client)[1])
  output_object_for_client <- relocate(output_object_for_client, rank)
  colnames(output_object_for_client)[colnames(output_object_for_client) == "X._TME_in_indication"] <- "%_TME_in_indication"
  
  
  head(output_object_for_client,nrow(output_object_for_client))
  
}


