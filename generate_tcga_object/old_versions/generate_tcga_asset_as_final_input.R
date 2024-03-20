# load env
library(ggpubr)
library(forcats)
library(dplyr)
library(boot, lib.loc = "/usr/lib/R/library")
devtools::load_all("~/capsule/code/cytoreason.ccm.pipeline/") 
BiocManager::install("hgu95av2.db")
n
library(hgu95av2.db)

# load the TCGA genes object (not final!!)
gene_exp_data <- get_outputs_dist("wf-1e2553df19")$output.rds
head(rownames(gene_exp_data),5) # the genes are located here as ENTREZID


### ~~ ### ~~ ### ~~ organize the object ~~ ### ~~ ### ~~ ###

# convert gene ENTREZID to SYMBOL
genes_as_symbols <- mapIds(org.Hs.eg.db, keys = rownames(gene_exp_data), keytype="ENTREZID", column = "SYMBOL")
head(genes_as_symbols, 10)
length(genes_as_symbols[is.na(genes_as_symbols)]) # 105 missing genes
gene_exp_data <- gene_exp_data[rownames(gene_exp_data) %in% names(genes_as_symbols[!is.na(genes_as_symbols)])]
genes_as_symbols <- genes_as_symbols[!is.na(genes_as_symbols)]
names(genes_as_symbols) <- NULL
rownames(gene_exp_data) <- genes_as_symbols
head(rownames(gene_exp_data),5) # the genes are located here as ENTREZID
rm(genes_as_symbols)

# concatenate fields to generate type_tme per sample
# The cluster annotation isn't the most updated - need to be changed later
gene_exp_data$type_tme <- paste0(gene_exp_data$tissue_io,"_",gene_exp_data$cluster_annotation)
table(gene_exp_data$type_tme)
length(unique(gene_exp_data$type_tme))
# [1] 30

# Add columns: n_per_indication, n_per_indication_TME, percent_tme
gene_exp_data <- data.frame(gene_exp_data)
gene_exp_data <- gene_exp_data %>%
  # generate n_per_indication and n_per_indication_TME columns
  dplyr::add_count(tissue_io, name = "n_per_indication") %>%
  dplyr::add_count(tissue_io, cluster_annotation, name = "n_per_indication_TME") %>%
  # generate percent_tme by the two variables above
  dplyr::mutate(percent_tme = (n_per_indication_TME / n_per_indication) * 100)


# Upload the final object to cyto-cc
tags <- list(list(name="owner", value="Ariel Simon"),
             list(name="data", value="indication prioritization input"),
             list(name="dataset", value="TCGA"))

wf <- run_function_dist(function(obj){return(obj)},
                        obj=gene_exp_data, tags = tags, image = 'eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:develop_latest') 
# Cyto-CC workflow: wf-3a3109e1f0