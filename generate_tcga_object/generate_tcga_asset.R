#### preparation ####
# libraries
library(cytoreason.cc.client)
library(checkmate)
library(plyr)
library(assertthat)
library(cytoreason.gx)
library(cytoreason.io)
library(cytoreason.shared.assets)
install.packages('cytoreason.ccm.pipeline')
library(cytoreason.ccm.pipeline)
# read the downloaded TCGA file
TCGA_all <- readRDS("~/capsule/code/TCGA_all.rds")
dim(TCGA_all)

# filter only for IO disease model samples
IO_tissues <- c('lung', 'prostate','bladder', 'breast', 'colon', 'head and neck')
io_tcga <- TCGA_all[,TCGA_all@phenoData@data$tissue %in% IO_tissues]

rm(TCGA_all); gc()

# correct the TCGA IO format for sample_id
io_tcga$sampleid <- gsub(".*_","",io_tcga$sampleid)
colnames(io_tcga) <- gsub(".*_","",colnames(io_tcga))
dim(io_tcga) # 4053 samples
names(io_tcga@phenoData@data)[2] <- "sample_id" # instead of "sampleid"
names(io_tcga@phenoData@data)[53] <- 'tissue_io'
rownames(io_tcga@phenoData@varMetadata)[2] <- "sample_id"
rownames(io_tcga@phenoData@varMetadata)[53] <- 'tissue_io'
# make sure that phenoData@data[2] is as the same size of phenoData@varMetadata[1] :)
dim(io_tcga@phenoData@data)
dim(io_tcga@phenoData@varMetadata)

# Upload the final object to cyto-cc
tags <- list(list(name="owner", value="Ariel Simon"),
             list(name="model", value="all 6 io models"),
             list(name="dataset", value="TCGA"))

wf <- run_function_dist(function(obj){return(obj)},
                        obj=io_tcga, tags = tags, image = 'eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:develop_latest') 
# Cyto-CC workflow: wf-1e2553df19