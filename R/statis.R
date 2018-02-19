#' STATIS, dual STATIS and SIR-STATIS for cubic data
#'
#' perform STATIS-like analyses of cubic data
#'
#' The package is well suited to explore time-course datasets, eventually in
#' relation with a target numeric variable (SIR approach).
#'
################################################################################
# main function
################################################################################
#' @title STATIS, dual STATIS and SIR-STATIS for cubic data.
#'
#' @description
#' \code{statis} performs the STATIS-like analyses of cubic data with the two
#' standard steps: inter-structure analysis (with the function
#' \code{\link{interStructure}}) and intra-structure analysis.
#'
#' @param data List; The list of tables to study.
#' @param method The method to apply ("sir", "classic" or "dual"). Default is
#' "sir".
#' @param y A numerical variable of interest (mandatory for "sir").
#' @param H Integer; The number of slices to cut y (mandatory for "sir").
#' @param scale Logical; Should the interstructure be computed with the cosine
#' (TRUE) or with the dot product (FALSE). Default is TRUE.
#' @return S3 object of class
#' \code{statisRes}: a list consisting of
#' \itemize{
#'    \item{parameters}{
#'    \itemize{
#'      \item{method}{ method used ("sir", "classic" or "dual")}
#'      \item{H}{ number of slices}
#'      \item{scale}{ logical indicating if the interstructure is analyzed with
#'      the cosine (TRUE) or the dot product (FALSE)}
#'      }
#'    }
#'    \item{svd}{
#'    \itemize{
#'      \item{values}{ eigenvalues of the generalized svd (square of singular
#'      values)}
#'      \item{P}{ left-singular vectors of the generalized svd}
#'      \item{Q}{ right-singular vectors of the generalized svd}
#'      }
#'    }
#'    \item{observations}{
#'    \itemize{
#'      \item{coord}{ if \code{method} is "sir" or "dual", matrix with the
#'      compromise coordinates of slices ("sir") or of observations (center of
#'      gravity of the table-specific coordinates). If \code{method} is
#'      "classic", matrix with the coordinates of the observations (factors).}
#'      \item{tcoord}{ if \code{method} is "sir" or "dual", list of matrices
#'      with the table-specific coordinates of slices ("sir") or of
#'      observations. If \code{method} is "classic", list of the coordinates of
#'      the partial (table-specific) factors}
#'      }
#'    }
#'    \item{variables}{
#'    \itemize{
#'      \item{cos2}{ if \code{method} is "sir" or "dual", matrix with the
#'      coordinates of the variables on the circle of correlation (also called
#'      squared cosinus). If \code{method} is "classic", matrix with the
#'      compromise coordinates of the variables on the circle of correlation
#'      (center of gravity of the table-specific coordinates)}
#'      \item{tcos2}{ if \code{method} is "sir" or "dual", list of matrices with
#'      the coordinates of the partial (table-specific) squared-cosinus. If
#'      \code{method} is "classic", list of the table-specific coordinates}
#'      }
#'    }
#'  }
#'  @export
statis <- function(data, method = c("sir", "classic","dual"), y = NULL,
                   H = NULL, scale = TRUE) {
  # perform standard checks
  method <- match.arg(method)
  all_nobs <- unlist(lapply(data, function(alist) nrow(alist)))
  all_nvars <- unlist(lapply(data, function(alist) nrow(alist)))
  if ((sum(all_nobs!=all_nobs[1]) > 0) || (sum(all_nvars!=all_nvars[1]) > 0))
    stop("the package requires cubic data (all data have the same dimension.")
  if (method == "sir") {
    if (is.null(y) || is.null(H))
      stop("arguments y and/or H are NULL (mandatory for 'sir').")
  }

  inter_struc <- interStructure(data, method, y, H, scale)
  transfo_data <- inter_struc$pdata$pdata

  # compute elements for the generalized SVD
  ## data
  if (method == "classic") {
    concatened_data <- do.call(cbind, transfo_data)
    concatened_data <- as.matrix(concatened_data)
  } else {
    concatened_data <- do.call(rbind, transfo_data)
    concatened_data <- as.matrix(concatened_data)
  }

  ## metric for the variables
  p <- ncol(transfo_data[[1]])
  if (method == "classic") {
    rootM <- diag(sqrt(inter_struc$tweights)) %x% diag(rep(1, p))
    invRM <- diag(1 / sqrt(inter_struc$tweights)) %x% diag(rep(1, p))
  } else {
    rootM <- diag(rep(1, p))
    invRM <- rootM
  }

  ## metric for the observations
  nobs <- nrow(transfo_data[[1]])
  if (method == "classic") {
    obs_weights <- rep(1, nobs)
    D <- diag(rep(1, nobs)) / nobs
    rootD <- diag(rep(1, nobs)) / sqrt(nobs)
    invRD <- diag(rep(1, nobs)) * sqrt(nobs)
  } else {
    obs_weights <- inter_struc$pdata$weights
    D <- diag(inter_struc$tweights) %x% diag(obs_weights) / nobs
    rootD <- diag(sqrt(inter_struc$tweights)) %x% diag(sqrt(obs_weights)) /
      sqrt(nobs)
    invRD <- diag(1 / sqrt(inter_struc$tweights)) %x%
      diag(1 / sqrt(obs_weights)) * sqrt(nobs)
  }

  # perform generalized SVD
  res_svd <- svd(rootD %*% concatened_data %*% rootM)
  P <- invRD %*% res_svd$u
  Q <- invRM %*% res_svd$v
  Sigma <- res_svd$d

  # coordinates
  ## representation of observations
  tsteps <- length(data)
  obs_coord <- P %*% diag(Sigma)
  if (method == "classic") {
    ### compromise
    coord <- obs_coord
    ### partial-factors
    tcoord <- list()
    for (ind in 1:tsteps) {
      sel_ind <- ((ind-1) * p + 1):(ind * p)
      tcoord[[ind]] <- concatened_data[ ,sel_ind] %*% Q[sel_ind, ]
    }
  } else {
    ### time-specific
    tcoord <- list()
    for (ind in 1:tsteps) {
      tcoord[[ind]] <- obs_coord[((ind-1) * nobs + 1):(ind * nobs), ]
    }
    ### compromise
    coord <- Reduce("+", mapply(function(tcoord, alpha) alpha * tcoord, tcoord,
                                inter_struc$tweights, SIMPLIFY = FALSE))
  }
  ### add names
  rownames(coord) <- rownames(transfo_data[[1]])
  tcoord <- lapply(tcoord, function(alist) {
    rownames(alist) <- rownames(transfo_data[[1]])
    alist
  })

  ## representation of variables
  var_coord <- t(concatened_data) %*% D %*% P
  st_dev <- sqrt(apply(sweep(concatened_data^2, 1, diag(D), "*"), 2, sum))
  var_coord <- sweep(var_coord, 1, st_dev, "/")
  if (method == "classic") {
    ### time-specific
    tcos2 <- list()
    for (ind in 1:tsteps) {
      sel_ind <- ((ind-1) * p + 1):(ind * p)
      tcos2[[ind]] <- var_coord[sel_ind, ]
    }
    ## compromise
    cos2 <- Reduce("+", mapply(function(tcoord, alpha) alpha * tcoord,
                               tcos2, inter_struc$tweights,
                               SIMPLIFY = FALSE))
  } else {
    ## compromise
    cos2 <- var_coord
    rownames(cos2) <- colnames(data[[1]])
    vcoord <- var_coord
    ## partial-loadings
    tcos2 <- list()
    for (ind in 1:tsteps) {
      sel_ind <- ((ind-1) * nobs + 1):(ind * nobs)
      tdata <- concatened_data
      tdata[-sel_ind, ] <- 0
      sd_var <- sqrt(apply(sweep(tdata^2, 1, diag(D), "*"), 2, sum))
      tcos2[[ind]] <- t(tdata) %*% D %*% P
      tcos2[[ind]] <- sweep(tcos2[[ind]], 1, sd_var, "/")
    }
  }
  ### add names
  rownames(cos2) <- colnames(data[[1]])
  tcos2 <- lapply(tcos2, function(alist) {
    rownames(alist) <- colnames(data[[1]])
    alist
  })

  # output
  parameters <- list("method" = method, "H" = H, "scale" = scale)

  res <- list("parameters" = parameters,
              "svd" = list("values" = res_svd$d^2, "P" = P, "Q" = Q),
              "observations" = list("coord" = coord, "tcoord" = tcoord),
              "variables" = list("cos2" = cos2, "tcos2" = tcos2),
              "interstructure" = inter_struc)

  class(res) <- "statisRes"
  return(res)
}

################################################################################
# Methods for objects of class statisRes
################################################################################
#' @S3method summary statisRes
summary.statisRes <- function(object,...) {
  cat("STATIS method with\n")
  if (object$parameters$scale) {
    scaling <- "scaled\n"
  } else scaling <- "not scaled\n"
  if (object$parameters$method == "sir") {
    cat("     ", object$parameters$method, ", ",
        object$parameters$H, "slices - ", scaling, "\n")
  } else
    cat("     ", object$parameters$method, " - ", scaling, "\n")
  cat("-----\n")
  print(object$interstructure)
  cat("-----\n")
  cat("     Compromise space\n")
  cat("Singular values\n", object$svd$values, "\n")
}

#' @S3method print statisRes
print.statisRes <- function(x,...) {
  summary(x)
  invisible(x)
}

