#' Confidence Interval for Multiresolution Functional ANOVA (MRFA) Model
#'
#' @description The function computes the confidence intervals of predicted responses (only works for linear regression model).
#'
#' @param object a class MRFA object estimated by \code{MRFA_fit}.
#' @param xnew a testing matrix with dimension \code{n_new} by \code{d} in which each row corresponds to a predictive location.
#' @param X input for \code{MRFA_fit}.
#' @param lambda a value. The default is \code{min(object$lambda)}.
#' @param conf.level a value specifying confidence level of the confidence interval. The default is 0.95.
#' @param parallel logical. If \code{TRUE}, apply function in parallel using parallel backend provided by foreach.
#' @param verbose logical. If \code{TRUE}, additional diagnostics are printed.
#'
#' @return
#' \item{lower bound}{a vector with length \code{n_new} displaying lower bound of predicted responses at locations \code{xnew}.}
#' \item{upper bound}{a vector with length \code{n_new} displaying upper bound of predicted responses at locations \code{xnew}.}
#' \item{conf.level}{as above.}
#'
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
#' }
#' @export confidence.MRFA
confidence.MRFA <- function(object, xnew, X, lambda = object$lambda, conf.level = 0.95, parallel = FALSE, verbose = FALSE){

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
  converge.tol          <-    object$converge.tol
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
  min.x                 <-    object$min.x
  scale.x               <-    object$scale.x

  X <- t((t(X) - min.x)/scale.x)    # scale X to [0,1]
  xnew <- t((t(xnew) - min.x)/scale.x)    # scale xnew to [0,1]

  n.xnew <- nrow(xnew)
  Phi <- basis_fun(active.group, X, gridpoint.ls, bandwidth.ls, parallel = parallel)
  Z <- cbind(1, Phi) ### active effects
  s <- ncol(Z)

  Q <- basis_fun(nonactive.group, X, gridpoint.ls, bandwidth.ls, parallel = parallel) ### nonactive effects
  Phi.train <- cbind(Z, Q)
  beta_hat <- cbind(beta_hat, matrix(0, nrow = nrow(beta_hat), ncol = ncol(Q)))

  Phi.test <- basis_fun(active.group, xnew, gridpoint.ls, bandwidth.ls, parallel = parallel)
  Phi.test <- cbind(1, Phi.test)
  Phi.test <- cbind(Phi.test, basis_fun(nonactive.group, xnew, gridpoint.ls, bandwidth.ls, parallel = parallel))
  p <- ncol(Phi.test)

  if(s < 30 & n > s){
    sigma <- sqrt(colSums((Y -  Phi.train %*% t(beta_hat))^2)/(n-s))
  }else{
    cv.out <- cv.MRFA(X = X, Y = Y, order = order, level = level, lambda = exp(seq(log(lambda + 1000), log(lambda), by = -0.01)),
                      converge.tol = converge.tol, nvar.max = s - 1, K = 5, plot.it = FALSE, parallel = parallel, verbose = verbose)
    sigma <- sqrt(cv.out$cv[length(cv.out$cv)])
  }

  if(parallel){
    if (foreach::getDoParWorkers() == 1) {
      # EXCLUDE COVERAGE START
      warning("No parallel backend registered", call. = TRUE)
      # EXCLUDE COVERAGE END
    }

    CI <- foreach::foreach(i = seq(n.xnew), .combine = rbind, .verbose = verbose, .packages = "glmnet") %dopar%{
      c1 <- Phi.test[i,, drop = FALSE]
      B <- rbind(c1, diag(1,length(c1))[2:(length(c1)),])
      #P <- Phi.train %*% solve(B)
      #Znew <- P[,1]
      #Qnew <- P[,2:ncol(P)]
      ### faster to compute the B inverse and P
      a <- rep(0, length(c1))
      a[1] <- 1/c1[1]
      a[2:length(a)] <- -c1[2:length(c1)]/c1[1]
      #B.inv <- rbind(a, diag(1,length(a))[2:(length(a)),])
      Znew <- a[1] * Phi.train[,1]
      Qnew <- Phi.train[,2:ncol(Phi.train)]
      for(j in 2:ncol(Phi.train)){
        Qnew[,j-1] <- a[j] * Phi.train[,1] + Phi.train[,j]
      }

      eta_hat <- B %*% t(beta_hat)

      x.lars <- crossprod(Qnew)
      y.lars <- crossprod(Qnew, Znew)

      m.out <- glmnet(x = x.lars, y = y.lars, intercept=FALSE, alpha = 1)
      cv.out <- cv.glmnet(x = x.lars, y = y.lars, lambda = seq(0, max(m.out$lambda), length = 100), intercept=FALSE, alpha = 1, nfolds = 5, parallel = parallel)
      cv.mode <- cv.out$lambda.min
      lasso.fit <- glmnet(x = x.lars, y = y.lars, lambda = cv.out$lambda, intercept=FALSE, alpha = 1)
      w <- coef(lasso.fit, s = cv.mode)[-1]

      U <- colSums((Y - Qnew %*% eta_hat[-1,]) * c(Znew - Qnew %*% w))
      b <- sum(Znew * (Znew - Qnew %*% w))
      c_alpha_1 <- U + sigma * qnorm(1-(1-conf.level)/2) * sqrt(b)
      c_alpha_2 <- U - sigma * qnorm(1-(1-conf.level)/2) * sqrt(b)

      if(b > 0){
        UCI <- c_alpha_1/b
        LCI <- c_alpha_2/b
      }else{
        UCI <- c_alpha_2/b
        LCI <- c_alpha_1/b
      }
      out <- rbind(UCI, LCI)
    }
    UCI <- CI[(1:n.xnew)*2-1,]
    LCI <- CI[(1:n.xnew)*2,]
  }else{

    UCI <- LCI <- matrix(0, nrow = n.xnew, ncol = 1)
    for(i in 1:n.xnew){
      c1 <- Phi.test[i,, drop = FALSE]
      B <- rbind(c1, diag(1,length(c1))[2:(length(c1)),])
      #P <- Phi.train %*% solve(B)
      #Znew <- P[,1]
      #Qnew <- P[,2:ncol(P)]
      ### faster to compute the B inverse and P
      a <- rep(0, length(c1))
      a[1] <- 1/c1[1]
      a[2:length(a)] <- -c1[2:length(c1)]/c1[1]
      #B.inv <- rbind(a, diag(1,length(a))[2:(length(a)),])
      Znew <- a[1] * Phi.train[,1]
      Qnew <- Phi.train[,2:ncol(Phi.train)]
      for(j in 2:ncol(Phi.train)){
        Qnew[,j-1] <- a[j] * Phi.train[,1] + Phi.train[,j]
      }

      eta_hat <- B %*% t(beta_hat)

      x.lars <- crossprod(Qnew)
      y.lars <- crossprod(Qnew, Znew)

      m.out <- glmnet(x = x.lars, y = y.lars, intercept=FALSE, alpha = 1)
      cv.out <- cv.glmnet(x = x.lars, y = y.lars, lambda = seq(0, max(m.out$lambda), length = 100), intercept=FALSE, alpha = 1, nfolds = 5, parallel = parallel)
      cv.mode <- cv.out$lambda.min
      lasso.fit <- glmnet(x = x.lars, y = y.lars, lambda = cv.out$lambda, intercept=FALSE, alpha = 1)
      w <- coef(lasso.fit, s = cv.mode)[-1]

      U <- colSums((Y - Qnew %*% eta_hat[-1,]) * c(Znew - Qnew %*% w))
      b <- sum(Znew * (Znew - Qnew %*% w))
      c_alpha_1 <- U + sigma * qnorm(1-(1-conf.level)/2) * sqrt(b)
      c_alpha_2 <- U - sigma * qnorm(1-(1-conf.level)/2) * sqrt(b)


      if(b > 0){
        UCI[i,] <- c_alpha_1/b
        LCI[i,] <- c_alpha_2/b
      }else{
        UCI[i,] <- c_alpha_2/b
        LCI[i,] <- c_alpha_1/b
      }
    }
  }

  return(list(UCI = c(UCI), LCI = c(LCI), conf.level = conf.level))
}
