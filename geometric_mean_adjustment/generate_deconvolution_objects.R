# load env
library(cytoreason.cc.client)
library(cytoreason.deconvolution)

# read the annotation data of our tool samples (3791 samples)
annotation_table <-  read_data(get_task_outputs("wf-2c26fb32e8", 0, 
                                                files_names_grepl_pattern = "tcga_annotation_data.csv")["tcga_annotation_data.csv"])

# read the raw expression matrix TPM values
TCGA_updated <- readRDS(get_workflow_outputs('wf-01d457dac0'))
expression_matrix <- assayDataElement(TCGA_updated, "abundance") # TPM values

# subset the entire expression_matrix (9164 samples) so only our 6 IO diseases samples will be there (3791 samples)
expression_matrix <- expression_matrix[, colnames(expression_matrix) %in% annotation_table$sample_id]
head(expression_matrix[1:5,1:5])
dim(expression_matrix)


# send the raw expression to deconvolution woth SO signature
deconve_raw_data <- run_function_dist(cytoreason.deconvolution::service_ct_contribution,
                                      data=ExpressionSet(expression_matrix),
                                      method = "cyto_marker",
                                      nperm = 0,
                                      signature_collection = "epithelial_tumor_myloid")
# Cyto-CC workflow: wf-b0067e79f7



### ~ Geometric mean normalization to the TPM counts ~ ###

# step 1 - to avoid negative or 0 values, add the minimal number in expression matrix + 0.01 to have only positive numbers
expression_matrix <- expression_matrix + abs(min(expression_matrix)) + 0.01
head(expression_matrix[1:5,1:5])

# step 2 - calc the geometric mean factor per sample
GZMA_PRF1 <- c('3001', '5551') # these entrez_ids are corresponding to the GZMA and PRF1 symbols
GZMA_PRF1_factor <- exp(colMeans(log(expression_matrix[GZMA_PRF1, ])))

# normalize each sample by the it's geometric mean factor 
expression_matrix <- sweep(expression_matrix, 2L, GZMA_PRF1_factor, FUN = "/")
head(expression_matrix[1:5,1:5])


# send the normalized expression to deconvolution woth SO signature
deconve_manipulated_data <- run_function_dist(cytoreason.deconvolution::service_ct_contribution,
                                              data=ExpressionSet(expression_matrix),
                                              method = "cyto_marker",
                                              nperm = 0,
                                              signature_collection = "epithelial_tumor_myloid")
# Cyto-CC workflow: wf-9d21bd3f11

