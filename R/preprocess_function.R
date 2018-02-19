################################################################################
# auxiliary functions
################################################################################

##### tools: compute the inverse square root of a matrix
#' @import MASS
#' @import expm
sqrtInvMatrix <- function (amatrix) {
  eig_matrix <- ginv(amatrix)
  return(sqrtm(eig_matrix))
}

#' @title Conditional expectation for the SIR method.
#'
#' @description
#' \code{condExpectation} computes the conditional expectation of X given y with
#' a slicing approach.
#'
#' @param data List; The list of tables to study.
#' @param y A numerical variable of interest.
#' @param H Integer; The number of slices to cut y.
#' @param scale Logicial; Should the data be scaled by the Frobenius norm of the
#' covariance matrix (TRUE) or not (FALSE)? Default is TRUE.
#' @return a list consisting of
#' \itemize{
#'  \item{pdata}{ the conditional expectations of X given y (list of lenght T)}
#'  \item{weights}{ the number of observations in every slides}
#'  }
#' @export
condExpectation <- function(data, y, H, scale = TRUE) {
  #centers data
  centered_data <- lapply(data, function(alist) {
    scale(alist, center = TRUE, scale = FALSE)
  })
  #computes the covmatrix for each table
  gamma <- lapply(data, function(alist) cov(alist))

  #computes gamma: covmatrix^(-1/2)
  sqrt_gamma <- lapply(gamma, function(alist) sqrtInvMatrix(alist))

  #computes Z..t: weights each t table by its gamma (after standardization)
  Z <- mapply(function(ttable, tgamma) as.matrix(ttable) %*% tgamma,
              centered_data, sqrt_gamma, SIMPLIFY = FALSE)

  slices <- cut(y, quantile(y, (0:H)/H), include.lowest = TRUE)
  weights <- table(slices)
  cond_exp <- lapply(Z, function(alist)
    apply(alist, 2, function(acol) tapply(acol, slices, mean, na.rm = TRUE))
  )
  cond_exp <- lapply(cond_exp, na.omit)

  if (scale) {
    D <- diag(weights / nrow(data[[1]]))
    scalar_prod <- lapply(cond_exp, function(alist)
      t(alist) %*% D %*% as.matrix(alist))
    frobenius_norm <- lapply(scalar_prod, function(alist)
      sqrt(sum(diag(alist %*% alist))))
    cond_exp <- mapply(function(value, fnorm) {
      value / sqrt(fnorm)
    }, cond_exp, frobenius_norm, SIMPLIFY = FALSE)
  }

  return (list("pdata" = cond_exp, "weights" = weights))
}
