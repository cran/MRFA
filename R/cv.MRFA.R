#' Compute K-fold cross-validated error for Multi-Resolution Functional ANOVA (MRFA) Model
#'
#' @description Computes the K-fold cross validated mean squared prediction error for multiresolution functional ANOVA model.
#'
#' @param X input for \code{MRFA_fit}.
#' @param Y input for \code{MRFA_fit}.
#' @param order input for \code{MRFA_fit}.
#' @param level input for \code{MRFA_fit}.
#' @param lambda lambda values at which CV curve should be computed.
#' @param K a positive integer specifying the number of folds.
#' @param plot.it logical. If \code{TRUE}, a CV curve will be shown. The default is \code{TRUE}.
#' @param parallel logical. If \code{TRUE}, apply cross-validation function in parallel using parallel backend provided by foreach. The default is \code{FALSE}.
#' @param verbose logical. If \code{TRUE}, additional diagnostics are printed. The default is \code{FALSE}.
#' @param ... additional arguments to \code{MRFA_fit}.
#'
#' @return
#' \item{lambda}{lambda values at which CV curve is computed.}
#' \item{cv}{the CV curve at each value of lambda.}
#' \item{cv.error}{the standard error of the CV curve}
#'
#' @seealso \code{\link{MRFA_fit}} for fitting a multiresolution functional ANOVA model.
#' @author Chih-Li Sung <iamdfchile@gmail.com>
#'
#' @import graphics
#' @examples \dontrun{
#'
#' #####             Testing function: GRAMACY & LEE (2009) function             #####
#' #####   Thanks to Sonja Surjanovic and Derek Bingham, Simon Fraser University #####
#' grlee09 <- function(xx)
#' {
#'   x1 <- xx[1]
#'   x2 <- xx[2]
#'   x3 <- xx[3]
#'   x4 <- xx[4]
#'   x5 <- xx[5]
#'   x6 <- xx[6]
#'
#'   term1 <- exp(sin((0.9*(x1+0.48))^10))
#'   term2 <- x2 * x3
#'   term3 <- x4
#'
#'   y <- term1 + term2 + term3
#'   return(y)
#' }
#'
#' library(MRFA)
#' #####   Training data and testing data   #####
#' set.seed(2)
#' n <- 100; n_rep <- 3; n_new <- 50; d <- 6
#' X.train <- matrix(runif(d*n), ncol = d)
#' X.train <- matrix(rep(X.train, each = n_rep), ncol = d)
#' Y.train <- apply(X.train, 1, grlee09)
#' Y.train <- Y.train + rnorm(n*n_rep, 0, 0.05)
#' X.test <- matrix(runif(d*n_new), ncol = d)
#' Y.test <- apply(X.test, 1, grlee09)
#'
#' #####   Fitting    #####
#' MRFA_model <- MRFA_fit(X.train, Y.train)
#'
#' #####   Computes the K-fold cross validated   #####
#' cv.out <- cv.MRFA(X.train, Y.train, K = 5, lambda = seq(0.01,3,0.1))
#'
#' #####   Prediction : CV  ######
#' lambda_cv <- cv.out$lambda[which.min(cv.out$cv)]
#' Y.pred <- predict(MRFA_model, X.test, lambda = lambda_cv)$y_hat
#' print(sqrt(mean((Y.test - Y.pred)^2)))
#' }
#' @export
#'

cv.MRFA <- function (X, Y, order = 10, level = 10, lambda = exp(seq(log(500), log(0.001), by = -0.01)),
                     K = 10, plot.it = TRUE, parallel = FALSE, verbose = FALSE, ...)
{
  ##### The code is modified from cv.lars (see the original package lars by Trevor Hastie)


  all.folds <- cv.folds(length(Y), K)
  lambda.min <- min(lambda)

  if(parallel){
    if (!requireNamespace("foreach", quietly = TRUE)) {
      # EXCLUDE COVERAGE START
      stop("foreach package required for parallel cv.MRFA operation",
           call. = FALSE)
      # EXCLUDE COVERAGE END
    }
    if (foreach::getDoParWorkers() == 1) {
      # EXCLUDE COVERAGE START
      warning("No parallel backend registered", call. = TRUE)
      # EXCLUDE COVERAGE END
    }
    residmat <- foreach::foreach(i = seq(K), .combine = cbind, .verbose = verbose) %dopar%{
      omit <- all.folds[[i]]
      fit <- MRFA_fit(X[-omit, , drop = FALSE], Y[-omit], order = order, level = level,
                      lambda.min = lambda.min, parallel = FALSE, verbose = verbose, ...)
      pred.object <- predict(fit, X[omit, , drop = FALSE], lambda = lambda, parallel = parallel)
      Y_hat <- pred.object$y_hat

      if (length(omit) == 1)
        Y_hat <- matrix(Y_hat, nrow = 1)

      out <- apply((Y[omit] - Y_hat)^2, 2, mean)
      pred.lambda <- sort(pred.object$lambda, decreasing = TRUE)
      if(which.min(pred.lambda) != length(pred.object$lambda)){
        out[lambda < min(pred.lambda)] <- NA
      }
      pred.lambda <- sort(pred.object$lambda)
      if(which.max(pred.lambda) != length(pred.object$lambda)){
        out[lambda > max(pred.lambda)] <- NA
      }

      return(out)
    }
  }else{
    residmat <- matrix(0, length(lambda), K)
    for (i in seq(K)) {
      omit <- all.folds[[i]]
      fit <- MRFA_fit(X[-omit, , drop = FALSE], Y[-omit], order = order, level = level,
                      lambda.min = lambda.min, parallel = parallel, verbose = verbose, ...)
      pred.object <- predict(fit, X[omit, , drop = FALSE], lambda = lambda, parallel = parallel)
      Y_hat <- pred.object$y_hat

      if (length(omit) == 1)
        Y_hat <- matrix(Y_hat, nrow = 1)

      residmat[, i] <- apply((Y[omit] - Y_hat)^2, 2, mean)
      pred.lambda <- sort(pred.object$lambda, decreasing = TRUE)
      if(which.min(pred.lambda) != length(pred.object$lambda)){
        residmat[lambda < min(pred.lambda), i] <- NA
      }
      pred.lambda <- sort(pred.object$lambda)
      if(which.max(pred.lambda) != length(pred.object$lambda)){
        residmat[lambda > max(pred.lambda), i] <- NA
      }

      if (verbose)
        cat("\n CV Fold", i, "\n\n")
    }
  }

  na.index <- apply(residmat, 1, function(x) any(is.na(x)))
  residmat <- residmat[!na.index,]
  lambda <- lambda[!na.index]
  cv <- apply(residmat, 1, mean)
  cv.error <- sqrt(apply(residmat, 1, var)/K)
  object <- list(lambda = lambda, cv = cv, cv.error = cv.error)
  if (plot.it){
    plotCVMRFA(object)
  }

  invisible(object)
}
