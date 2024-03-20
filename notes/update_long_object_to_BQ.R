
library(SingleCellExperiment)
library(cytoreason.io)
library(cytoreason.ccm.pipeline)
library('cytoreason')
library(pkgmaker)
library(cytoreason.datasource)
library(naturalsort)
# library(cytoreason.single.cell)

library(cytoreason)
library(devtools)

library(cytoreason.cc.client)

# upload to BQ
wf_id='wf-bc9435edd8'
dataset_name = 'p00_io_tcga' 
IMAGE <- 'eu.gcr.io/cytoreason/cd-py-bigquery:develop_latest'


command <- paste("python /app/cytobigquery/exec_service.py indication_prioritization", wf_id, dataset_name, "--verbose -sc -ti=eu.gcr.io/cytoreason/cd-py-bigquery:develop_latest -tm=124Mi")
res <-
  run_command_dist(
    command,
    image = IMAGE,
    memory_request = "25Gi",
    force_execution = FALSE,
    replace_image_tags = TRUE
  )

# https://cyto-cc.cytoreason.com/workflow/wf-987256b92e