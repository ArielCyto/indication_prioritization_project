#' ---
#' title: indication prioritization -- ccm-run
#' author: Model Team
#' fontsize: 9pt
#' output: html_document
#' params:
#'   signature_collection: "pan_cancer_v1"
#'   ccm_pipeline_version: "master@0.57.0" #"master@0.57.2"
#'   configuration_file: "ccm-metadata.csv"
#'   tags:
#'     disease_model: "io models"
#'     tissue: "breast"
#'     condition: "io models"
#'     project: "indication prioritization"
#' ---
#'
#' ----------------
#
#' # Setup {.tabset}
#' ## Packages
#+ setup, include = FALSE

knitr::opts_chunk$set(comment = "#>", collapse = TRUE)
install.packages('unixtools', repos = 'http://www.rforge.net/')
unixtools::set.tempdir("/scratch")


# install.packages("cytoreason.integration")
library("cytoreason.integration")
# install.packages("cytoreason.deconvolution")
library("cytoreason.deconvolution")
install.packages("cytoreason.ccm.pipeline")
library(cytoreason.ccm.pipeline)
# loads ccm pipeline package
library(checkmate)
library(assertthat)
#devtools::load_all("~/capsule/code/cytoreason.ccm.pipeline/") 
library(cytoreason.io)
library(cytoreason.project)
# prepare output directory
library(cytoreason.shared.assets)
library(cytoreason.datasource)
setwd("~/capsule/code/CRC_disease_model/notes/indication_prioritization/")

IMAGE_v3 = "eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:SCAL_128_integrate_so_classifier_to_ccm_0.57.14"
# setwd("~/capsule/code/analysis-disease-models/notes/lung_cancer/lung_cancer")
OUTDIR <- cytoreason.shared.assets::get_output_dir(SCRIPT <- "generate-tcga-ccm.R")

#' ## Parameters
#+ params, echo = FALSE
params <- cytoreason.project::get_script_parameters(file = SCRIPT, use_script_defaults = TRUE)
str(params)
# extract parameters
PIPELINE_VERSION <- getObjectElement(params, "ccm_pipeline_version")
assert_string(ccm_cyto_cc_image_name(image = PIPELINE_VERSION))
#
TAGS <- getObjectElement(params, "tags", class = "list")
assert_list(TAGS, types = "character", names = "unique")
assert_subset(c("disease_model", "tissue", "condition"), names(TAGS))
#
SIGNATURE_COLLECTION <- getObjectElement(params, "signature_collection")
try(SIGNATURE_COLLECTION <- AssetData(SIGNATURE_COLLECTION), silent = TRUE)
SIGNATURE_COLLECTION


CONFIG_FILE = read.csv("ccm-metadata.csv")
CONFIG_FILE <- getObjectElement(params, "configuration_file")
# this checks the config syntax and sets the asset_id so that it pulls from the validator (it needs to exist)

# Here you might need to restart R session if you get:
# Error in (function (classes, fdef, mtable) :
#             unable to find an inherited method for function 'save_data' for signature '"data.frame", "character"'
#install.packages("cytoreason.ccm.pipeline")
#library(cytoreason.ccm.pipeline)
ccm_prepare_configuration(CONFIG_FILE, overwrite = TRUE)
library(cytoreason.shared.assets)
file.copy(CONFIG_FILE, file.path(OUTDIR, paste0("model-metadata.", tools::file_ext(CONFIG_FILE))), overwrite = TRUE)
DT::datatable(read_data(CONFIG_FILE))
OBJECTS <- ccm_stage_prepare_dataset_collection(CONFIG_FILE)

IMAGE <- IMAGE_v3
TAGS$comment <- "second trial indication prioritization"
(wf <- ccm_api_generate_disease_model(CONFIG_FILE, #qc = list(qc_measures=c("adjustment", "covariates_outlier", "duplications", "meta-analysis",
                                                    #                       qc_memory_request="50Gi")),
                                      prepare_data = list(memory_request = "2Gi"),
                                      dataset = list(
                                        cell_contribution = list(signature_collection = SIGNATURE_COLLECTION)
                                      ),
                                      adjustment_models = list("1" = list(c(1),c(2), c(1, 2))),# c("memory B cell", "CD4-positive, alpha-beta T cell", "epithelial cell", "fibroblast", "endothelial cell", "Macro_LYVE1", "fat cell", "Macro_IL1B"))), ##c(1, 2), c(1, 3))),
                                      model = .skip("cell_specific_differences"),
                                      n_pc = 10,
                                      # qc = T,
                                      image = IMAGE,
                                      tags = TAGS, memory_request = "50Gi"))
# save_data(wf, file.path(OUTDIR, "ccm_workflow.qs"))



####### Step 2 - Generate text output #######

library(cytoreason.cc.client)
library(cytoreason.ccm.pipeline)

disease_name <- 'p00_io_breast_cancer'
ccm_wf_id <-"wf-ede58a8618"  #"wf-7aeb37b59b" 
memory_request <- '32Gi'
memory_request_downstream <- '[20Gi]'
classifier <- "V3"
special_details <- "check subtypes TCGA - with IO flag"
TAGS$comment <- "check subtypes TCGA - with IO flag"
wf <- ccm_api_save_data(AssetData(ccm_wf_id),
                        DirectoryTXT(), #bio_qc = F,
                        image = ccm_cyto_cc("eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:develop_latest",
                                            save_data.api = memory_request_downstream),
                        # image = ccm_cyto_cc("eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:feature_ENGINE_2044_support_json_config_0.57.17",
                        #                     save_data.api = memory_request_downstream),
                        # 
                        memory_request=memory_request,
                        tags=list(list(name="disease", value=disease_name),
                                  list(name="ccm-wf-id", value=ccm_wf_id),
                                  list(name="classifier", value=classifier),
                                  list(name="comment", value=TAGS$comment),
                                  #  list(name="project", value=TAGS$project),
                                  list(name="special_details", value=special_details)
                        )
)

# save_data -- Tue Oct 24 13:04:06 2023: wf-27e4bc3bf1


####### Step 3 - upload text output #######

library(cytoreason.cc.client)

#########
# image #
#########
task_image = 'eu.gcr.io/cytoreason/cd-py-bigquery:develop_latest'
memory_request = '64Mi'

#####################
# ARGS (positional) #
#####################
service = 'ccm'
dataset = "p00_io_breast_cancer" ## !!!!!!!! Pay attention to the DB name !!!!!!!! #
wf_id = "wf-51d3f9b1b8" #"wf-2b915f046f" # "wf-27e4bc3bf1"


###################
# ARGS (optional) #
###################
verbose = '--verbose'
ti = sprintf('-ti=%s', task_image)
tm = sprintf('-tm=%s', memory_request)


###########
# command #
###########
command <- sprintf('python /app/cytobigquery/exec_service.py %s %s %s %s %s %s',
                   service, wf_id, dataset, verbose, ti, tm)
############
# workflow #
############
tags=list(list(name="de-process", value='bigquery-upload'),
          list(name="service", value=service),
          list(name="export_workflow", value=wf_id),
          list(name="target_dataset", value=dataset)#,
          # list(name="project", value=TAGS$project)
)
task_env_vars = list(list("name"="DE_PROCESS","value"="BigQuery_upload")
)
wf <- run_command_dist(command,
                       outdir = "./output/",
                       image = task_image ,
                       task_env_vars = task_env_vars,
                       tags=tags,
                       memory_request = memory_request)
wf
# Cyto-CC workflow: wf-21c54815ed
# Cyto-CC workflow: wf-93e8b7a3a0
