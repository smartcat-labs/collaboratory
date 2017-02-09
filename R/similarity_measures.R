# ============================================================
# Similarity measures for sparse matrices.
# ============================================================

#' Calculates correlations between columns of two sparse matrices.
#'
#' @param X (dgCMatrix)
#' @param Y (dgCMatrix)
#' @returns Matrix of correlations.
cal_cor <- function(X, Y){

  availX <- X!=0
  availY <- Y!=0
  
  # TODO: to optimize further
  X<- t(as(t(X) - colMeans(X), "dgCMatrix"))
  Y<- t(as(t(Y) - colMeans(Y), "dgCMatrix"))

  R <- crossprod(X,Y)
  N <- crossprod(X^2, availY)
  M <- crossprod(availX, Y^2)

  cor <- R
  cor@x <- cor@x/((N@x^0.5) * (M@x^0.5))

  cor
}

#' Calculates cosine between columns of two sparse matrices.
#'
#' @param X (dgCMatrix)
#' @param Y (dgCMatrix)
#' @returns Matrix of cosine measures.
cal_cos <- function(X, Y){
  
  ones <- rep(1,nrow(X))		
  means <- drop(crossprod(X^2, ones)) ^ 0.5
  diagonal <- Diagonal( x = means^-1 )
  X <- X %*% diagonal
  
  ones <- rep(1,nrow(Y))		
  means <- drop(crossprod(Y^2, ones)) ^ 0.5
  diagonal <- Diagonal( x = means^-1 )
  Y <- Y %*% diagonal
  
  crossprod(X, Y)
}
