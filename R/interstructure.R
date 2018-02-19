################################################################################
# inter-structure analysis
################################################################################

#' @title Inter-structure analysis.
#'
#' @description
#' \code{interStructure} perfoms the inter-structure analysis and returns its
#' result (table weights, compromise similarity matrix between tables,...). The
#' function is also called by \code{\link{statis}}. You should only use it when
#' wanting to perform only the inter-structure analysis.
#'
#' @param data List; The list of tables to study.
#' @param method The method to apply ("sir", "classic" or "dual"). Default is
#' "sir".
#' @param y A numerical variable of interest (mandatory for "sir").
#' @param H Integer; The number of slices to cut y (mandatory for "sir").
#' @param scale Logical; Should the interstructure be computed with the cosine
#' (TRUE) or with the dot product (FALSE). Default is TRUE.
#' @return S3 \code{interstructure} object; a list consisting of
#' \itemize{
#'  \item{objects}{ representative objects according to the chosen method}
#'  \item{corr}{ matrix used for the eigendecomposition}
#'  \item{values}{ eigenvalues of the decomposition of the RV coefficient
#'  or correlation coefficient matrix}
#'  \item{vectors}{ eigenvectors of the decomposition of the RV coefficients
#'  or correlation coefficient matrix}
#'  \item{consensus}{ interstructure consensual matrix}
#'  \item{tweights}{ computed weights for each table}
#'  \item{pdata}{ a list with two elements: \code{pdata} the preprocessed data
#'  (list with length T) used to compute the scalar products (centered and
#'  eventually scaled; conditional expectations for "sir" as produced by
#'  \link{condExpectation}) and \code{weights} the weights associated with every
#'   observation.}
#'  }
#' @export
interStructure <- function(data, method = c("sir", "classic", "dual"), y = NULL,
                           H = NULL, scale = TRUE) {
  method <- match.arg(method)
  all_nobs <- unlist(lapply(data, function(alist) nrow(alist)))
  all_nvars <- unlist(lapply(data, function(alist) nrow(alist)))
  if ((sum(all_nobs!=all_nobs[1]) > 0) || (sum(all_nvars!=all_nvars[1]) > 0))
    stop("the package requires cubic data (all data have the same dimension).")

  if (method == "sir") {
    if (is.null(y) || is.null(H))
      stop("arguments y and/or H are NULL (mandatory for 'sir').")

    #computes the conditional expectation matrix for each ttable
    cond_matrix <- condExpectation(data, y, H, scale)
    #computes frequencies for each H class of y for each t table
    D <- diag(cond_matrix$weights / nrow(data[[1]]))
    #computes the distance between variables on the transformed data
    scalar_prod <- lapply(cond_matrix$pdata, function(alist)
      t(alist) %*% D %*% as.matrix(alist))
    scalar_prod <- lapply(scalar_prod, function(x) {
      colnames(x) <- rownames(x) <- colnames(data[[1]])
      return (x)})
  } else if (method %in% c("classic", "dual")) {
    data <- lapply(data, scale, center = TRUE, scale = FALSE)
    if (method == "dual") {
      scalar_prod <- lapply(data, function(alist)
        t(alist) %*% as.matrix(alist) / nrow(alist))
      scalar_prod <- lapply(scalar_prod, function(x) {
        colnames(x) <- rownames(x) <- colnames(data[[1]])
        return (x)})
    } else {
      scalar_prod <- lapply(data, function(alist) as.matrix(alist) %*% t(alist))
      scalar_prod <- lapply(scalar_prod, function(x) {
        colnames(x) <- rownames(x) <- rownames(data[[1]])
        return (x)})
    }
    if (scale) {
      frobenius_norm <- lapply(scalar_prod, function(alist)
        sqrt(sum(diag(alist %*% alist))))
      data <- mapply(function(tdata, fnorm) {tdata / sqrt(fnorm)}, data,
                     frobenius_norm, SIMPLIFY = FALSE)
      scalar_prod <- mapply(function(tdata, fnorm) {tdata / fnorm}, scalar_prod,
                            frobenius_norm, SIMPLIFY = FALSE)
    }
  }

  #computes the RV matrix
  timeS <- length(scalar_prod)
  cor_time <- outer(1:timeS, 1:timeS, FUN = Vectorize(function(i,j) {
    sum(diag(as.matrix(scalar_prod[[i]]) %*% scalar_prod[[j]]))
  }))

  #PCA
  resPCA <- eigen(cor_time)
  resPCA$vectors[ ,1] <- abs(resPCA$vectors[ ,1])
  consensus <- Reduce("+", mapply(function(mat, weights) {
    mat * weights / sum(resPCA$vectors[ ,1])
  }, scalar_prod, resPCA$vectors[ ,1], SIMPLIFY = FALSE))

  res <- list("objects" = scalar_prod, "corr" = cor_time,
              "values" = resPCA$values, "vectors" = resPCA$vectors,
              "tweights" = resPCA$vectors[,1] / sum(resPCA$vectors[ ,1]),
              "consensus" = consensus)

  #add (preprocessed) original data to the output
  if (method == "sir") {
    res$pdata <- cond_matrix
  } else if (method %in% c("classic", "dual")) {
    #note: for dual, this part is only valid when all data have the same number
    #of observations which is not mandatory in the general case
    res$pdata <- list("pdata" = data, "weights" = rep(1, nrow(data[[1]])))
  }

  class(res) <- "interstructure"
  return(res)
}

################################################################################
# Methods for objects of class interstructure
################################################################################
#' @S3method summary interstructure
summary.interstructure <- function(object,...) {
  cat("     Interstructure\n")
  cat("Compromise quality\n",
      round(object$values[1] * 100 / sum(object$values), 2), "%\n")
  cat("Eigenvalues\n", object$values, "\n")
  cat("Weights\n", object$tweights, "\n")
}

#' @S3method print interstructure
print.interstructure <- function(x,...) {
  summary(x)
  invisible(x)
}
