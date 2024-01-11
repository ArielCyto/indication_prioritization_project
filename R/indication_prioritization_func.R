# Goal: to check that the gene or genes who were provided are available in our data
# Input: gene_list - gene or genes in SYMBOL format. examples: one gene 'CD4'. list of genes c('BRCA1','BRCA2', 'HOXB13')
# Output: submit - logical variable that determine if gene_list is valid (and than the analysis starts automatically). If no.
# the user get a message and the analysis does not start.

check_input_genes <- function(gene_list){
  # check if the genes are available:
  submit = FALSE
  for (gene in gene_list){
    gene_status <- gene %in% colnames(expression_matrix)
    if (gene_status){
      submit = TRUE
    }else{
      print(paste0("gene ",gene," is an invalid input, please remove it and try again" ))
      submit = FALSE
      optional_alternative_name <- colnames(expression_matrix)[(grepl(toupper(gene),colnames(expression_matrix)))]
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
  bo <- by(data_table, data_table$group, function(i) boot(i, gene_stat, R = repeats))
  bo_ci <- t(sapply(bo, function(x) {
    bo_ci <- boot.ci(x, type = "basic")
    # bo_ci keeps 3 fields per type_tme: t0, low_ci and high_ci
    c(t0 = bo_ci$t0, CI_low = bo_ci$basic[4], CI_high = bo_ci$basic[5])}))
  bo_ci <- as.data.frame(bo_ci)
  bo_ci$group <- rownames(bo_ci)
  
  
  data_table <- left_join(data_table, bo_ci, by = "group")
  return(data_table)
}

generate_output_object <- function(data_table_after_boots){
  
  
  output_object <-
    data_table_after_boots %>%
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
  return(output_object)
}

print_ranking_plot <- function(output_object,gene_list, rank_by){
  output_object$group = with(output_object, reorder(group, output_object[[rank_by]]))
  output_object$rank_by <- output_object[[rank_by]] # syntax for using ggplot
  print(output_object %>%
          slice_head(n = 10) %>%
          ggplot(aes(x = rank_by, y = group)) + ggtitle(paste0("Ranking by ", paste(gene_list, collapse = ', ')))+
          labs(x = as.character(rank_by))+
          geom_point(alpha = 0.8, size = 3, aes(colour = factor(TME))) +
          geom_errorbar(aes(xmin = CI_low, xmax = CI_high),
                        width = .2, size = .8,
                        position = "identity", colour = "darkgray"))
}


print_ranking_table <- function(output_object, rank_by = rank_by){
  output_object_for_client <- data.frame(output_object[order(output_object[[rank_by]], decreasing = TRUE),])
  output_object_for_client$rank <- 1:(dim(output_object_for_client)[1])
  output_object_for_client <- relocate(output_object_for_client, rank)
  colnames(output_object_for_client)[colnames(output_object_for_client) == "X._TME_in_indication"] <- "%_TME_in_indication"
  head(output_object_for_client,nrow(output_object_for_client))
}

# main function
rank_groups_by_genes <- function(gene_list, repeats = 5, rank_by = 'CI_low'){ # need to add option for specific TME/ indication
 
  # check if the input list is valid
  submit = check_input_genes(gene_list)  

  # run analysis
  if (submit == TRUE){
    
    # step 1 - per sample (row), calculate the selected gene_list average expression
    if (length(gene_list) > 1)
    {expression_matrix$gene_list_average <- rowMeans(expression_matrix[,gene_list])
    }else{
      expression_matrix$gene_list_average <- expression_matrix[,gene_list]
    }
    
    # work with smaller object of gene_exp_data with all the necessary columns
    data_table <- annotation_table
    data_table$gene_list_average <- expression_matrix$gene_list_average
    
    
    ## BOOTSTRAPING - calculate the CI per type_tme per collection (based on the effector cells average)
    set.seed(42)
    data_table <- get_statistics(data_table,repeats = repeats)

    # generate output object
    output_object <- generate_output_object(data_table)
    
    # visualize results
    # print results plot
    print_ranking_plot(output_object, gene_list, rank_by)
    
    # print results numeric
    print_ranking_table(output_object, rank_by)
    
    
  }
}
