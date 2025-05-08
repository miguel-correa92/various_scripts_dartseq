# Calculate Modified Roger's Distance from allele frequency data
mrd_parallel <- function(mat, n_cores = NULL) {
  require(magrittr)
  require(foreach)
  require(doParallel)
  require(dplyr)
  
  if (!is.matrix(mat)) {
    stop('input not a matrix!')
  }
  
  n_cols <- ncol(mat)
  mrdmat <- matrix(NA, nrow = n_cols, ncol = n_cols)
  
  # Set up parallel backend
  if (is.null(n_cores)) {n_cores <- parallel::detectCores() - 5}
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  paths <- .libPaths()
  
  results <- foreach(samp = 1:(n_cols - 1), .combine = rbind, .packages = 'magrittr' ) %dopar% {

    if (samp != n_cols) {
      temp <- colMeans(
        (mat[, samp] - mat[, samp:n_cols])^2 +
          ((1 - mat[, samp]) - (1 - mat[, samp:n_cols]))^2,
        na.rm = TRUE
      ) %>% sqrt * 1 / sqrt(2)
      temp <- c(rep(NA, samp - 1), temp)
      temp
    } else {
      rep(NA, n_cols)
    } 
  }
  
  # Stop parallel backend
  stopCluster(cl)
  
  # Fill the results into the matrix
  for (samp in 1:(n_cols - 1)) {
    mrdmat[samp, ] <- results[samp, ]
  }
  
  colnames(mrdmat) <- colnames(mat)
  rownames(mrdmat) <- colnames(mat)
  diag(mrdmat) <- 0
  
  return(mrdmat)
}
