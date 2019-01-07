#' Prediction of Multi-Resolution Functional ANOVA (MRFA) Model
#'
#' @description The function computes the predicted responses.
#'
#' @param object a class MRFA object estimated by \code{MRFA_fit}.
#' @param xnew a testing matrix with dimension \code{n_new} by \code{d} in which each row corresponds to a predictive location.
#' @param lambda a value, or vector of values, indexing the path. The default is \code{object$lambda}.
#' @param parallel logical. If \code{TRUE}, apply function in parallel in \code{ldply} using parallel backend provided by foreach.
#' @param ... for compatibility with generic method \code{predict}.
#'
#' @return
#' \item{lambda}{as above.}
#' \item{coefficients}{coefficients with respect to the basis function value.}
#' \item{y_hat}{a matrix with dimension \code{n_new} by \code{length(lambda)} displaying predicted responses at locations \code{xnew}.}
#'
#' @seealso \code{\link{MRFA_fit}} for fitting a multiresolution functional ANOVA model.
#' @author Chih-Li Sung <iamdfchile@gmail.com>
#'
#' @examples \dontrun{
#'
#' #####             Testing function: OTL circuit function                     #####
#' #####  Thanks to Sonja Surjanovic and Derek Bingham, Simon Fraser University #####
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
#' library(MRFA)
#' #####   Training data and testing data   #####
#' set.seed(2)
#' n <- 1000; n_new <- 100; d <- 6
#' X.train <- matrix(runif(d*n), ncol = d)
#' Y.train <- apply(X.train, 1, otlcircuit)
#' X.test <- matrix(runif(d*n_new), ncol = d)
#' Y.test <- apply(X.test, 1, otlcircuit)
#'
#' #####   Fitting    #####
#' MRFA_model <- MRFA_fit(X.train, Y.train, verbose = TRUE)
#'
#' #####   Prediction   ######
#' Y.pred <- predict(MRFA_model, X.test, lambda = min(MRFA_model$lambda))$y_hat
#' print(sqrt(mean((Y.test - Y.pred)^2)))
#' }
#' @export predict.MRFA
#' @export

predict.MRFA <- function(object, xnew, lambda = object$lambda, parallel = FALSE, ...){

  ##### modified from predict.lars (see original package by Trevor Hastie)

  if (is.MRFA(object) == FALSE) {
    stop("The object in question is not of class \"MRFA\" \n")
  }
  if (is.matrix(xnew) == FALSE) {
    xnew <- as.matrix(xnew)
  }

  coefficients          <-    object$coefficients
  active.group          <-    object$active.group
  gridpoint.ls          <-    object$gridpoint.ls
  bandwidth.ls          <-    object$bandwidth.ls
  index                 <-    object$index
  active_group.index    <-    object$active_group.index
  min.x                 <-    object$min.x
  scale.x               <-    object$scale.x
  k                     <-    object$k
  standardize.d         <-    object$standardize.d

  ###################       Setting       ###################
  if(standardize.d) xnew <- t((t(xnew) - min.x)/scale.x)    # scale xnew to [0,1]
  p <- ncol(coefficients)
  n.test <- nrow(xnew)
  Phi <- basis_fun(active.group, xnew, gridpoint.ls, bandwidth.ls, k, parallel)
  Phi <- cbind(rep(1,nrow(Phi)), Phi)

  ##################       Prediction      ###################
  row.index <- is.element(index, active_group.index) | is.na(index)
  betas <- t(object$coefficients[row.index, ])
  dimnames(betas) <- list(NULL, dimnames(betas)[[2]])
  kp <- dim(betas)
  kk <- kp[1]
  p <- kp[2]

  lambdas <- object$lambda
  lambda[lambda > max(lambdas)] <- max(lambdas)
  lambda[lambda < min(lambdas)] <- min(lambdas)
  sbeta <- lambdas

  sfrac <- (lambda - sbeta[1])/(sbeta[kk] - sbeta[1])
  sbeta <- (sbeta - sbeta[1])/(sbeta[kk] - sbeta[1])
  usbeta <- unique(sbeta)
  useq <- match(usbeta, sbeta)
  sbeta <- sbeta[useq]
  betas <- betas[useq, , drop = FALSE]
  coord <- approx(sbeta, seq(sbeta), sfrac)$y
  left <- floor(coord)
  right <- ceiling(coord)
  newbetas <- ((sbeta[right] - sfrac) * betas[left, , drop = FALSE] +
                 (sfrac - sbeta[left]) * betas[right, , drop = FALSE])/(sbeta[right] -
                                                                          sbeta[left])
  newbetas[left == right, ] <- betas[left[left == right], ]
  y_hat <- drop(Phi %*% t(newbetas))

  return(list(lambda = lambda, coefficients = newbetas, y_hat = y_hat))
}
