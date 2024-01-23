service_ssgsea <- function(data, collection, alpha = 0.25, 
                           scale = TRUE, normalize = TRUE,
                           drop = TRUE) {
  # extract expression data if not directly a matrix
  if( !is.matrix(data) ) data <- exprs(data)
  
  assert_matrix(data, row.names = "unique")
  assert_list(collection, types = "character", names = "unique")
  assert_number(alpha, lower = 0, upper = 1)
  assert_flag(scale)
  assert_flag(normalize)
  assert_flag(drop)
  
  # Ranks for genes
  message("ssGSEA - data (%s x %s | %s gene sets) ", 
           nrow(data), ncol(data), 
           length(collection),
           appendLF = FALSE)
  fids <- rownames(data)
  n_features <- nrow(data)
  
  # powered ranks ----
  message("> ranks ", appendLF = FALSE)
  R_alpha <- matrixStats::colRanks(data, preserveShape = T, ties.method = 'average') ^ alpha
  # incidence ----
  message("> incidence (sparse) ", appendLF = FALSE)
  # NOTE: we use a sparse format which speeds-up ~ 7x the computations
  j <- lapply(lapply(collection, match, fids), function(x) x[!is.na(x)])
  i <- unlist(mapply(rep.int, seq_along(j), lengths(j), SIMPLIFY = FALSE))
  gs_I <- sparseMatrix(i, unlist(j), 
                       dims = c(length(collection), length(fids)),
                       dimnames = list(names(collection), fids)) # -> class ngCMatrix
  gs_size <- rowSums(gs_I)
  assert_true(identical(rownames(gs_I), names(collection)))
  # cumsum factors ----
  message("> cumsum factors ", appendLF = FALSE)
  idx <- seq(n_features)
  cumsum_f <- n_features + 1 - apply(R_alpha, 2L, function(r) match(idx, order(r, decreasing = TRUE)))
  # negative part ----
  message("> negative part ", appendLF = FALSE)
  # NOTE: we compute negative part as the complementary of the sum of the 
  # factors within the gene set because this calculation is better making 
  # use of the sparseness of gs_I
  neg_mat <- sweep(-(gs_I %*% cumsum_f), 2L, colSums(cumsum_f), "+")
  neg_mat <- sweep(neg_mat, 1L, n_features - gs_size, "/")
  # positive part ----
  message("> positive part ", appendLF = FALSE)
  pos_mat <- (gs_I %*% (R_alpha * cumsum_f)) / (gs_I %*% R_alpha)
  es <- pos_mat - neg_mat
  
  # cast to dense matrix
  es <- as.matrix(es)
  # scale --- 
  # Normalize by gene number
  if (scale){
    message("> scale ", appendLF = FALSE)
    es <- es / n_features
    
  }
  
  # normalize ----
  # Normalize by absolute diff between max and min
  if (normalize){
    message("> normalize ", appendLF = FALSE)
    es = es / diff(range(es, na.rm = TRUE))
    
  }
  message("[OK]")
  
  # fix values for empty gene-sets: set to NA
  if( !drop ) es[gs_size == 0, ] <- NA_real_
  else es <- es[gs_size > 0, , drop = FALSE]
  
  # Prepare output
  colnames(es) <- colnames(data)
  
  return(es)
  
  
}