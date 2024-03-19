# load env
library(cytoreason.cc.client)
library(cytoreason.deconvolution)

# read the data
annotation_table <-  read_data(get_task_outputs("wf-2c26fb32e8", 0, 
                                                files_names_grepl_pattern = "tcga_annotation_data.csv")["tcga_annotation_data.csv"])
TCGA_updated <- readRDS(get_workflow_outputs('wf-01d457dac0'))
expression_matrix <- assayDataElement(TCGA_updated, "abundance") # TPM values
head(expression_matrix[1:5,1:5])

# subset expression_matrix so only our 6 IO diseases samples will be there
expression_matrix <- expression_matrix[, colnames(expression_matrix) %in% annotation_table$sample_id]
dim(expression_matrix)
head(expression_matrix[1:5,1:5])
expression_matrix['3001',1:5]

# to avoid negative or 0 values, we will add the minimal number in the table + 0.01 to have only positive numbers
expression_matrix <- expression_matrix + abs(min(expression_matrix)) + 0.01
# head(expression_matrix[1:5,1:5])

# calc the geometric mean per sample
GZMA_PRF1 <- c('3001', '5551')
GZMA_PRF1_scaling <- exp(colMeans(log(expression_matrix[GZMA_PRF1, ])))

# normalize all the matrix by the geometric mean of the sample 
expression_matrix <- sweep(expression_matrix, 2L, GZMA_PRF1_scaling, FUN = "/")
head(expression_matrix[1:5,1:5])
expression_matrix['3001',1:5]



# toy_data <- expression_matrix[, 1:5] 
# deconve_manipulated_toy_data <- run_function_dist(cytoreason.deconvolution::service_ct_contribution,
#                                                   data=ExpressionSet(toy_data),
#                                                   method = "cyto_marker",
#                                                   nperm = 0,
#                                                   signature_collection = "epithelial_tumor_myloid")
# 
# 
# # Cyto-CC workflow: wf-11da0ac6ec

deconve_manipulated_data <- run_function_dist(cytoreason.deconvolution::service_ct_contribution,
                                                  data=ExpressionSet(expression_matrix),
                                                  method = "cyto_marker",
                                                  nperm = 0,
                                                  signature_collection = "epithelial_tumor_myloid")
# Cyto-CC workflow: wf-9d21bd3f11


# send also the original expression matrix to deconve
expression_matrix_raw <- assayDataElement(TCGA_updated, "abundance") # TPM values
# subset expression_matrix so only our 6 IO diseases samples will be there
expression_matrix_raw <- expression_matrix_raw[, colnames(expression_matrix_raw) %in% annotation_table$sample_id]
dim(expression_matrix_raw)

deconve_raw_data <- run_function_dist(cytoreason.deconvolution::service_ct_contribution,
                                              data=ExpressionSet(expression_matrix_raw),
                                              method = "cyto_marker",
                                              nperm = 0,
                                              signature_collection = "epithelial_tumor_myloid")
# Cyto-CC workflow: wf-b0067e79f7
