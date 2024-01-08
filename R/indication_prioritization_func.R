# Goal: to check that the gene or genes who were provided are available in our data
# Input: gene_list - gene or genes in SYMBOL format. examples: one gene 'CD4'. list of genes c('BRCA1','BRCA2', 'HOXB13')
# Output: submit - logical variable that determine if gene_list is valid (and than the analysis starts automatically). If no.
# the user get a message and the analysis does not start.

check_input_genes <- function(gene_list){
  # check if the genes are available:
  submit = FALSE
  for (gene in gene_list){
    gene_status <- gene %in% colnames(gene_exp_data)
    if (gene_status){
      submit = TRUE
    }else{
      print(paste0("gene ",gene," is an invalid input, please remove it and try again" ))
      submit = FALSE
      optional_alternative_name <- colnames(gene_exp_data)[(grepl(toupper(gene),colnames(gene_exp_data)))]
      if (length(optional_alternative_name) > 0){
        print("optional alternatives: ")
        print(optional_alternative_name)
      }

      break
    }
  }
  return(submit)  
}


# Goal: define the statistic method which will be used as input for the bootstrapping.
# This function is a permanent object for get_statistics() function.
gene_stat <- function(d, i) median(d[i, "gene_list_average"])


get_statistics <- function(data_table, repeats = repeats) {
  bo <- by(data_table, data_table$type_tme, function(i) boot(i, gene_stat, R = repeats))
  bo_ci <- t(sapply(bo, function(x) {
    bo_ci <- boot.ci(x, type = "basic")
    # bo_ci keeps 3 fields per type_tme: t0, low_ci and high_ci
    c(t0 = bo_ci$t0, CI_low = bo_ci$basic[4], CI_high = bo_ci$basic[5])}))
  bo_ci <- as.data.frame(bo_ci)
  bo_ci$type_tme <- rownames(bo_ci)
  
  
  data_table <- left_join(data_table, bo_ci, by = "type_tme")
  return(data_table)
}

generate_output_object <- function(data_table_after_boots){

  output_object <- 
    data_table_after_boots %>%
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
  
  return(output_object)
}

print_ranking_plot <- function(output_object,gene_list){
  output_object$type_tme = with(output_object, reorder(type_tme, score))
  print(output_object %>%
          slice_head(n = 10) %>%
          ggplot(aes(x = score, y = type_tme)) + ggtitle(paste0("Ranking by ", paste(gene_list, collapse = ', ')))+
          geom_point(alpha = 0.8, size = 3, aes(colour = factor(cluster_annotation))) +
          geom_errorbar(aes(xmin = CI_low, xmax = CI_high),
                        width = .2, size = .8,
                        position = "identity", colour = "darkgray"))
}


print_ranking_table <- function(output_object, top = 30){
  output_object_for_client <- data.frame(output_object[order(output_object$score, decreasing = TRUE),])
  head(output_object_for_client,top)
}





# main function
rank_groups_by_genes <- function(gene_list, repeats = 5){ # need to add option for spesific TME/ indication
 
  # check if the input list is valid
  submit = check_input_genes(gene_list)  

  # run analysis
  if (submit == TRUE){
    
    # step 1 - per sample (row), calculate the selected gene_list average expression
    if (length(gene_list) > 1)
    {gene_exp_data$gene_list_average <- rowMeans(gene_exp_data[,gene_list])
    }else{
      gene_exp_data$gene_list_average <- gene_exp_data[,gene_list]
    }
    
    # work with smaller object of gene_exp_data with all the necessary columns
    data_table <- data.frame(type_tme = gene_exp_data$type_tme, cluster_annotation = gene_exp_data$cluster_annotation,
                             tissue_io = gene_exp_data$tissue_io,
                             gene_list_average = gene_exp_data$gene_list_average,
                             n_per_indication_TME = gene_exp_data$n_per_indication_TME,
                             percent_tme = gene_exp_data$percent_tme)
    
    
    
    ## BOOTSTRAPING - calculate the CI per type_tme per collection (based on the effector cells average)
    set.seed(42)
    data_table <- get_statistics(data_table,repeats = repeats)

    # generate output object
    output_object <- generate_output_object(data_table)
    
    # visualize results
    # print results plot
    print_ranking_plot(output_object, gene_list)
    
    # print results numeric
    print_ranking_table(output_object,30)
    
    
  }
}
