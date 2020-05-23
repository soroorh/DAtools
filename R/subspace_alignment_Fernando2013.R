#' Subspace alignment
#'
#' @export

# need two function for here: one for determining d_max and one for the actual alignmnet


# determine d_max ------
# PCAs <- prcomp(t(CPM), scale. = TRUE, center = TRUE)
# PCAt <- prcomp(t(CPM_PDX), scale. = TRUE, center = TRUE)
#
# ns <- ncol(CPM)
# nt <- ncol(CPM_PDX)
# n_min <- min(ns, nt)
#
# s_eigenval_diffs <- diff(PCAs$sdev^2,1)
# t_eigenval_diffs <- diff(PCAt$sdev^2,1)
#
# d_lim <- min(length(s_eigenval_diffs), length(t_eigenval_diffs))
#
# seq_lambda_diffs  <- Map("min", abs(t_eigenval_diffs), abs(s_eigenval_diffs)[1:d_lim])
#
# sigma <- 0.1
# gamma <- 10^5
# B <- 1
# lower_bound <- (1+sqrt(log(2/sigma)/2))*((16*((1:d_lim)^(1.5))*B)/(gamma*sqrt(n_min)))
#
# (unlist(seq_lambda_diffs) >= lower_bound) # maximum dimension that satisfies the inequality is 11 == d_lim
#
# d_max <- 11
#
#
# # subspace alignment  ------
#
# # d <- 6
# Xs <- PCAs$rotation[,1:d]
# Xt <- PCAt$rotation[,1:d]
# Xa <- Xs%*%t(Xs)%*%Xt
# Sa <- t(CPM)%*%Xa
# Tt <- t(CPM_PDX)%*%Xt
