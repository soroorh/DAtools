#' Subspace alignment
#' @param  sourceFeatures  source dataset
#' @param  targetFeatures  target dataset
#' @param  d               dimension of the subspace
#' @param  sigma           related to theoratical lower bound. 0.1 by default
#' @param  gamma           related to theoratical lower bound. 10^5 by default
#' @param  B               related to theoratical lower bound. 1 by default
#'
#'
#' @return list object containing transformed (aligned) source and target data, d, lower bound, and differences in eigen values
#' @export
alignSubspace <- function(sourceFeatures, targetFeatures, d=NULL, sigma=0.1, gamma=10^5,
                          B=1){

  # need two function for here: one for determining d_max and one for the actual alignmnet
  # determine d_max ------
  PCAs <- prcomp(sourceFeatures, scale. = TRUE, center = TRUE)
  PCAt <- prcomp(targetFeatures, scale. = TRUE, center = TRUE)

  ns <- nrow(sourceFeatures)
  nt <- nrow(targetFeatures)
  n_min <- min(ns, nt)

  s_eigenval_diffs <- diff(PCAs$sdev^2,1)
  t_eigenval_diffs <- diff(PCAt$sdev^2,1)

  d_lim <- min(length(s_eigenval_diffs), length(t_eigenval_diffs))

  seq_lambda_diffs  <- Map("min", abs(t_eigenval_diffs), abs(s_eigenval_diffs)[1:d_lim])

  lower_bound <- (1+sqrt(log(2/sigma)/2))*((16*((1:d_lim)^(1.5))*B)/(gamma*sqrt(n_min)))

  #(unlist(seq_lambda_diffs) >= lower_bound) # maximum dimension that satisfies the inequality is 11 == d_lim

  d_max <- max(which(unlist(seq_lambda_diffs) >= lower_bound))


  # subspace alignment  ------
  if(is.null(d)){
    d <- d_max
    cat("Using optimal value of d=",d," that satisfies the lower bound to determine the subspace \n")
  }
  Xs <- PCAs$rotation[,1:d]
  Xt <- PCAt$rotation[,1:d]
  Xa <- Xs%*%t(Xs)%*%Xt
  Sa <- sourceFeatures%*%Xa
  Tt <- targetFeatures%*%Xt

  return(list(alignedSource = Sa, alignedTarget = Tt,
              d=d,
              lowerBound = lower_bound,
              eigenvalDiffs = seq_lambda_diffs))

}
