#' CORAL
#' @param sourceFeatures source data
#' @param targetFeatures target data
#'
#' @export
coral <- function(sourceFeatures, targetFeatures){
  cat("computing covariance of source data \n")
  C_s <- cov(sourceFeatures) + diag(dim(sourceFeatures)[2])
  cat("computing covariance of target data \n")
  C_t <- cov(targetFeatures) + diag(dim(targetFeatures)[2])

  cat("computing whitenings \n")
  # whitening source
  D_s <- sourceFeatures %*% pracma::sqrtm(C_s)$Binv

  # re-coloring with target covariance
  D_s_adjusted <- D_s %*% pracma::sqrtm(C_t)$Binv

  cat("done \n")
  # output adjusted source data
  return(D_s_adjusted)
}
