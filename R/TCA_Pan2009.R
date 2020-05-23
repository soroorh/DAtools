#' TCA
#'
#' @param sourceFeatures source data
#' @param targetFeatures target data
#' @param m number of components
#' @param sigma width of RBF kernel
#' @param mu trade-off parameter
#' @export

# To do:
## 1.In bioinformatics n << P = number of genes, hence the number of latent spaces is limied by the number of samples.
## 2.need to code this properly. also consider returning median of all distances to guide optimal values of sigma
## related to 2.,Yair et al (2019) [Optimal transport on manifold of SPD matrices] says sigma^2 in RBF kernel is typically
## taken as the median of all eucleadian distances.

tca <- function(sourceFeatures, targetFeatures,m, sigma=1, mu = 0.1 ){
  # sourceFeatures: numeric matrix of source data. Samples in rows, features in columns
  # targetFeatures: numeric matrix of target data. Samples in rows, features in columns
  # m: number of transfer components (number of leading eigenvalues)
  # sigma: parameter sigma for gaussian kernel
  # mu: trade-off parameter
  if (m >= dim(sourceFeatures)[2]) stop('m is larger than the total number of features')

  # compute K
  Kss <- rdetools::rbfkernel(X=sourceFeatures, sigma = sigma, Y=sourceFeatures)
  Kst <- rdetools::rbfkernel(X=sourceFeatures, sigma = sigma, Y=targetFeatures)
  Kts <- rdetools::rbfkernel(X=targetFeatures, sigma = sigma, Y=sourceFeatures)
  Ktt <- rdetools::rbfkernel(X=targetFeatures, sigma = sigma, Y=targetFeatures)
  K <- rbind(cbind(Kss, Kst), cbind(Kts,Ktt))

  # compute L
  n1 <- dim(sourceFeatures)[1]
  n2 <- dim(targetFeatures)[1]
  L <- matrix(0, n1+n2, n1+n2)
  L[1:n1,1:n1] <- 1/(n1^2)
  L[(n1+1):(n1+n2),(n1+1):(n1+n2)] <- 1/(n2^2)
  L[L==0] <- -1/(n1*n2)

  # compute mu and H
  H <- diag(n1+n2) - matrix(1, n1+n2, n1+n2)/(n1+n2)

  # calculate m leading eigenvalues and corresponding eigenvectors
  I <- diag(n1+n2)
  M <- solve(I + mu*K%*%L%*%K)%*%(K%*%H%*%K)
  s <- svd(M)
  V <- s$v
  D <- s$d
  W <- V[,order(D, decreasing = TRUE)[1:m]]
  colnames(W) <- paste0("TC", 1:ncol(W))
  rownames(W) <- c(rownames(sourceFeatures), rownames(targetFeatures))

  # calculate the resulting kernel matrix
  Kr <- K%*%W%*%t(W)%*%K

  # compute recommended sigma squared (sigma2)
  all_dists <- pdist::pdist(X=sourceFeatures, Y=targetFeatures)
  opt.sigma2 = median(as.matrix(all_dists))

  return(list(W=W, K = Kr, opt.sigma2 = round(opt.sigma2,0)))
  }
