#' Confidence Interval for Multiresolution Functional ANOVA (MRFA) Model
#'
#' @description The function computes the confidence intervals of predicted responses (only works for linear regression model).
#'
#' @param object a class MRFA object estimated by \code{MRFA_fit}.
#' @param xnew a testing matrix with dimension \code{n_new} by \code{d} in which each row corresponds to a predictive location.
#' @param X input for \code{MRFA_fit}.
#' @param lambda a value. The default is \code{min(object$lambda)}.
#' @param var.estimation a character string specifying the estimation method for variance. "rss" specifies residual sum of squares, "cv" specifies a cross-validation method with \code{K} fold, and "posthoc" specifies a post-hoc estimation method. The default is "rss".
#' @param w.estimation a character string specifying the estimation method for weights w. "cv" specifies a cross-validation method with \code{K} fold, and "nugget" specifies a least square error method with nugget=\code{nugget}. The default is "cv".
#' @param K a positive integer specifying the number of folds.
#' @param nugget a value specifying the nugget value for \code{w.estimation}. The default is 1e-6. It only works when \code{w.estimation="nugget"}.
#' @param conf.level a value specifying confidence level of the confidence interval. The default is 0.95.
#' @param parallel logical. If \code{TRUE}, apply function in parallel using parallel backend provided by foreach.
#' @param verbose logical. If \code{TRUE}, additional diagnostics are printed.
#'
#' @return
#' \item{lower bound}{a vector with length \code{n_new} displaying lower bound of predicted responses at locations \code{xnew}.}
#' \item{upper bound}{a vector with length \code{n_new} displaying upper bound of predicted responses at locations \code{xnew}.}
#' \item{conf.level}{as above.}
#'
#' @details When The details about \code{var.estimation} and \code{w.estimation} can be seen in Sung et al. (2017+).
#' @seealso \code{\link{MRFA_fit}} for fitting of a multi-resolution functional ANOVA model; \code{\link{predict.MRFA}} for prediction of a multi-resolution functional ANOVA model.
#' @author Chih-Li Sung <iamdfchile@gmail.com>
#'
#' @import glmnet
#' @examples \dontrun{
#'
#' #####             Testing function: OTL circuit function                      #####
#' #####   Thanks to Sonja Surjanovic and Derek Bingham, Simon Fraser University #####
#' otlcircuit <- function(xx)
#' {
#'   Rb1  <- 50   + xx[1] * 100
#'   Rb2  <- 25   + xx[2] * 45
#'   Rf   <- 0.5  + xx[3] * 2.5
#'   Rc1  <- 1.2  + xx[4] * 1.3
#'   Rc2  <- 0.25 + xx[5] * 0.95
#'   beta <- 50   + xx[6] * 250
#'
#'   Vb1 <- 12*Rb2 / (Rb1+Rb2)
#'   term1a <- (Vb1+0.74) * beta * (Rc2+9)
#'   term1b <- beta*(Rc2+9) + Rf
#'   term1 <- term1a / term1b
#'
#'   term2a <- 11.35 * Rf
#'   term2b <- beta*(Rc2+9) + Rf
#'   term2 <- term2a / term2b
#'
#'   term3a <- 0.74 * Rf * beta * (Rc2+9)
#'   term3b <- (beta*(Rc2+9)+Rf) * Rc1
#'   term3 <- term3a / term3b
#'
#'   Vm <- term1 + term2 + term3
#'   return(Vm)
#' }
#'
#'
#'
#' library(MRFA)
#' #####   training data and testing data   #############
#' set.seed(2)
#' n <- 100; n_new <- 10; d <- 6
#' X.train <- matrix(runif(d*n), ncol = d)
#' Y.train <- apply(X.train, 1, otlcircuit)
#' X.test <- matrix(runif(d*n_new), ncol = d)
#' Y.test <- apply(X.test, 1, otlcircuit)
#'
#' #####   Fitting    #####
#' MRFA_model <- MRFA_fit(X.train, Y.train)
#'
#' #####   Prediction   ######
#' Y.pred <- predict(MRFA_model, X.test, lambda = min(MRFA_model$lambda))$y_hat
#' print(sqrt(mean((Y.test - Y.pred)^2)))
#'
#' ### confidence interval ###
#' conf.interval <- confidence.MRFA(MRFA_model, X.test, X.train, lambda = min(MRFA_model$lambda))
#' print(conf.interval)
#' }
#' @export confidence.MRFA
confidence.MRFA <- function(object, xnew, X, lambda = object$lambda, conf.level = 0.95,
                            var.estimation = c("rss", "cv", "posthoc")[1], w.estimation = c("cv", "nugget")[1], K = 5, nugget = 1e-6,
                            parallel = FALSE, verbose = FALSE){

  if (is.MRFA(object) == FALSE) {
    stop("The object in question is not of class \"MRFA\" \n")
  }
  if (is.matrix(xnew) == FALSE) {
    xnew <- as.matrix(xnew)
  }

  if (length(lambda) > 1) {
    stop("length(lambda) > 1 is not allowed.\n")
  }

  if(object$model@name != "Linear Regression Model"){
    stop("The function only works for linear regression model.\n")
  }

  Y                     <-    object$y
  order                 <-    object$order
  level                 <-    object$level
  lambda.min            <-    object$lambda.min
  converge.end          <-    object$converge.end
  nvar.max              <-    object$nvar.max
  coefficients          <-    object$coefficients
  active.group          <-    object$active.group
  candidate.group       <-    object$candidate.group
  gridpoint.ls          <-    object$gridpoint.ls
  bandwidth.ls          <-    object$bandwidth.ls
  index                 <-    object$index
  active_group.index    <-    object$active_group.index
  n                     <-    nrow(X)
  nonactive.group       <-    candidate.group[-active_group.index]
  beta_hat              <-    predict(object, xnew, lambda, parallel)$coefficients
  nvar.end              <-    length(beta_hat)
  min.x                 <-    object$min.x
  scale.x               <-    object$scale.x

  X <- t((t(X) - min.x)/scale.x)          # scale X to [0,1]
  xnew <- t((t(xnew) - min.x)/scale.x)    # scale xnew to [0,1]
  n.xnew <- nrow(xnew)

  unique.ls <- unique.matrix(active.group, beta_hat, gridpoint.ls)
  unique.beta <- unique.ls$unique.beta
  unique.active.group <- unique.ls$unique.active.group

  Phi <- basis_fun(unique.active.group, X, gridpoint.ls, bandwidth.ls, parallel = parallel)  ### active effects
  s <- ncol(Phi) + 1

  unique.nonactive.group <- lapply(nonactive.group, function(x) x[length(x)])
  Q <- basis_fun(unique.nonactive.group, X, gridpoint.ls, bandwidth.ls, parallel = parallel) ### nonactive effects
  Phi.train <- cbind(Phi, Q)
  beta_hat <- cbind(matrix(unique.beta, nrow = 1), matrix(0, nrow = nrow(beta_hat), ncol = ncol(Q)))

  Phi.test <- basis_fun(unique.active.group, xnew, gridpoint.ls, bandwidth.ls, parallel = parallel)
  Phi.test <- cbind(Phi.test, basis_fun(unique.nonactive.group, xnew, gridpoint.ls, bandwidth.ls, parallel = parallel))
  p <- ncol(Phi.test) + 1

  if(s >= n & var.estimation == "rss"){
    warning("var.estimation is forced to be cv since s >= n")
    var.estimation <- "cv"
  }
  if(var.estimation == "rss"){
    sigma <- sqrt(colSums((Y -  Phi.train %*% t(beta_hat[,-1, drop=FALSE]) - beta_hat[1,1])^2)/(n-s))
  }else if (var.estimation == "cv"){
    cv.out <- cv.MRFA(X = X, Y = Y, order = order, level = level, lambda = exp(seq(log(lambda + 1000), log(lambda), by = -0.01)),
                      converge.tol = converge.end, nvar.max = nvar.end - 1, K = K, plot.it = FALSE, parallel = parallel, verbose = verbose)
    sigma <- min(sqrt(cv.out$cv))
  }else if(var.estimation == "posthoc"){
    sigma <- post_sigma(Phi.train, Y, beta_hat, conf.level = 0.95,
                        var.estimation, w.estimation, K, nugget, parallel, verbose)
  }else{
    stop("var.estimation setting is incorrect.")
  }

  XtX <- t(Phi.train) %*% Phi.train
  Xsum <- matrix(apply(Phi.train, 2, sum), nrow = 1)
  res <- Y - Phi.train %*% t(beta_hat[,-1, drop=FALSE])

  if(parallel){
    if (foreach::getDoParWorkers() == 1) {
      # EXCLUDE COVERAGE START
      warning("No parallel backend registered", call. = TRUE)
      # EXCLUDE COVERAGE END
    }

    CI <- foreach::foreach(i = seq(n.xnew), .combine = rbind, .verbose = verbose, .packages = "glmnet") %dopar%{

      c1 <- cbind(1, Phi.test[i,, drop = FALSE])
      B <- rbind(c1, diag(1,length(c1))[2:(length(c1)),])
      a <- rep(0, length(c1))
      a[1] <- 1/c1[1]
      a[2:length(a)] <- -c1[2:length(c1)]/c1[1]
      Znew <- a[1] * matrix(1, ncol = 1, nrow = n)

      xXsum <- t(c1[,-1, drop = FALSE]) %*% Xsum
      x.lars <- XtX - (xXsum + t(xXsum)) / c1[1] + n * t(c1[,-1, drop = FALSE]) %*% c1[,-1, drop = FALSE] / c1[1]^2
      y.lars <- t(Xsum - n * c1[,-1, drop = FALSE] / c1[1])

      if(w.estimation == "cv"){
        m.out <- glmnet(x = x.lars, y = y.lars, intercept = FALSE, alpha = 1)
        cv.out <- cv.glmnet(x = x.lars, y = y.lars, lambda = exp(seq(log(1e-6), log(max(m.out$lambda)), length = 100)), intercept=FALSE, alpha = 1, nfolds = K, parallel = parallel)
        cv.mode <- cv.out$lambda.min
        w <- coef(cv.out$glmnet.fit, s = cv.mode)[-1]
      }else{
        w <- solve(x.lars + diag(1e-6, ncol(x.lars)), y.lars)
      }

      res2 <- Znew - Phi.train %*% w
      U <- colSums((res - sum(a[-1] * beta_hat[,-1])) * c(res2  - sum(a[-1] * w)))
      b <- sum(res2 - sum(a[-1] * w))

      c_alpha_1 <- U + sigma * qnorm(1-(1-conf.level)/2) * sqrt(b)
      c_alpha_2 <- U - sigma * qnorm(1-(1-conf.level)/2) * sqrt(b)

      UCI <- c_alpha_1/b
      LCI <- c_alpha_2/b

      out <- rbind(UCI, LCI)
    }
    UCI <- CI[(1:n.xnew)*2-1,]
    LCI <- CI[(1:n.xnew)*2,]
  }else{
    UCI <- LCI <- matrix(0, nrow = n.xnew, ncol = 1)
    for(i in 1:n.xnew){
      c1 <- cbind(1, Phi.test[i,, drop = FALSE])
      B <- rbind(c1, diag(1,length(c1))[2:(length(c1)),])
      a <- rep(0, length(c1))
      a[1] <- 1/c1[1]
      a[2:length(a)] <- -c1[2:length(c1)]/c1[1]
      Znew <- a[1] * matrix(1, ncol = 1, nrow = n)

      xXsum <- t(c1[,-1, drop = FALSE]) %*% Xsum
      x.lars <- XtX - (xXsum + t(xXsum)) / c1[1] + n * t(c1[,-1, drop = FALSE]) %*% c1[,-1, drop = FALSE] / c1[1]^2
      y.lars <- t(Xsum - n * c1[,-1, drop = FALSE] / c1[1])

      if(w.estimation == "cv"){
        m.out <- glmnet(x = x.lars, y = y.lars, intercept = FALSE, alpha = 1)
        cv.out <- cv.glmnet(x = x.lars, y = y.lars, lambda = exp(seq(log(1e-6), log(max(m.out$lambda)), length = 100)), intercept=FALSE, alpha = 1, nfolds = K, parallel = parallel)
        cv.mode <- cv.out$lambda.min
        w <- coef(cv.out$glmnet.fit, s = cv.mode)[-1]
      }else{
        w <- c(solve(x.lars + diag(nugget, ncol(x.lars)), y.lars))
      }

      res2 <- Znew - Phi.train %*% w
      U <- colSums((res - sum(a[-1] * beta_hat[,-1])) * c(res2  - sum(a[-1] * w)))
      b <- sum(res2 - sum(a[-1] * w))

      c_alpha_1 <- U + sigma * qnorm(1-(1-conf.level)/2) * sqrt(b)
      c_alpha_2 <- U - sigma * qnorm(1-(1-conf.level)/2) * sqrt(b)

      UCI[i,] <- c_alpha_1/b
      LCI[i,] <- c_alpha_2/b
    }
  }

  return(list(UCI = c(UCI), LCI = c(LCI), conf.level = conf.level))
}
