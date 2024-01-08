### ~ Indication prioritization TCGA object - version 1 ~ ###
# TCGA asset (6 IO indications) : wf-1e2553df19
# CCM : wf-2f3002f651
# Genes tables : wf-218b83ea42

# env
library(cytoreason.cc.client)
library(checkmate)
library(plyr)
library(assertthat)
library(cytoreason.gx)
library(cytoreason.io)
library(cytoreason.shared.assets)
devtools::load_all("~/capsule/code/cytoreason.ccm.pipeline/")
library(plyr)

# # read TCGA CCM object and extract relevant the gene tables
# ccm_fit <- as_ccm_fit(read_asset("wf-2f3002f651"))
# gx_diff_res <- read_asset(ccm_fit$datasets$TCGA$model$io_cancer_TME__TCGA__io_cancer@analysis$gx_diff$fits$bulk)
# gx_diff_res <- gx_diff_res$stats$io_cancer_TME
# gx_diff_res <- gx_diff_res[1:29] # we want only the relevant terms of 29 groups - prostate_cancer:inflamed_mild_stromal is not found
# 
# # Upload the final tables to Cyto-cc
# tags <- list(list(name="owner", value="Ariel Simon"),
#              list(name="project", value="indication prioritization"),
#              list(name="data", value="TCGA - IO"))
# wf <- run_function_dist(function(obj){return(obj)},
#                         obj=gx_diff_res, tags = tags, image = "eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:develop_latest") 
# # Cyto-CC workflow: wf-218b83ea42

#########################################################################

genes_table <- get_outputs_dist('wf-218b83ea42')$output.rds
  
# define genes ranking function (estimation based)
get_genes_ranking <- function(term_list, gene){
  score_df <- list()
  for (term_name in term_list){
    subtable <- genes_table[term_name]
    subtable <- subtable[[1]]
    score_df[[term_name]] <- subtable[subtable$feature_id == gene,]$estimate
  }
  score_df <- data.frame(do.call("rbind",score_df))
  colnames(score_df) <- "gene_score"
  score_df$group <- rownames(score_df)
  rownames(score_df) <- NULL
  
  score_df <- score_df[order(score_df$gene_score, decreasing = TRUE),]
  return(score_df)
}


### Use the function above for scoring the groups by different genes
terms_list <- as.character(names(genes_table))

# print top results for cxcl8 gene (3576)
rank_3576 <- get_genes_ranking(terms_list,3576)
head(rank_3576,10)

# print top results for CD4 gene (920) - immune marker
rank_920 <- get_genes_ranking(terms_list,920)
head(rank_920,29)


# print top results for ACTA2 gene (59) - stromal marker
rank_59 <- get_genes_ranking(terms_list,59)
head(rank_59,7)

# print top results for ACTA2 gene (59) - stromal marker, within lung cancer specific
rank_59_lung <- get_genes_ranking(terms_list[grepl(pattern = 'lung',x = terms_list)],59)
head(rank_59_lung)
                  