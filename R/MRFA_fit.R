#' Fit a Multi-Resolution Functional ANOVA (MRFA) Model
#'
#' @description The function performs the multi-resolution functional ANOVA (MRFA) approach.
#'
#' @param X a design matrix with dimension \code{n} by \code{d}.
#' @param Y a response vector of size \code{n}.
#' @param weights a vector of observation weights.
#' @param order a positive integer specifying the highest order of interactions that can be entertained in the model. The default is 10.
#' @param level a positive integer specifying the highest resolution level that can be entertained in the model. The default is 10.
#' @param lambda.min a positive value specifying the minimum penalty value to be performed before the convergence criterion is met.
#' @param converge.tol convergence tolerance. It converges when relative difference with respect to function value (penalized likelihood) is smaller than the tolerance. The default is 1e-10.
#' @param nvar.max maximum number of non-zero variables.
#' @param pen.norm a character string specifying the type of penalty norm for group lasso to be computed. "2" or 2 specifies 2-norm, and "N" specifies native norm. The default is "2".
#' @param model an object of class specifying other models. \code{LinReg()} (default) fits a linear regression, \code{LogReg()} fits a logistic regression, and \code{PoissReg()} fits a Poisson regression.
#' @param parallel logical. If \code{TRUE}, apply function in parallel in \code{ldply} using parallel backend provided by foreach.
#' @param verbose logical. If \code{TRUE}, additional diagnostics are printed.
#'
#' @details A multi-resolution functional ANOVA (MRFA) model targets a low resolution representation of a low order functional ANOVA, with respect to strong effect heredity, to form an accurate emulator in a large-scale and high dimensional problem. This function fits an MRFA model using a modified group lasso algrithm. One can consider the loss function \deqn{\frac{1}{n}\sum_{i=1}^n\left(y_i-\sum_{|u|=1}^{D_{\rm max}}\sum_{r=1}^{R_{\rm max}}\sum_{k=1}^{n_u(r)}\beta_u^{rk}\varphi_u^{rk}(x_{iu})\right)^2+\lambda\sum_{|u|=1}^{D_{\rm max}}\sum_{r=1}^{R_{\rm max}}\sqrt{N_u(r)\sum_{v\subseteq u}\sum_{s\le r}\sum_{k=1}^{n_v(s)}(\beta_v^{sk})^2},} where \eqn{\varphi_u^{rk}(x_{iu})} is the basis function with resolution level \eqn{r} and with dimension \eqn{u\subset\{1,2,\ldots,d\}}, and \eqn{D_{\rm max}} and \eqn{R_{\rm max}} respectively are the maximal orders of functional ANOVA and multi-resolution level, which are indicated by \code{order} and \code{level}.
#'
#' The group lasso path along the penalty parameter \eqn{\lambda} is given by the function, where the \eqn{\lambda_{\max}} is automatically given and \eqn{\lambda_{\min}} is given by users, which is indicated by \code{lambda.min}. The group lasso algrithm is implemented via the modifications to the source code of the \code{grplasso} package (Meier, 2015).
#'
#' \code{lambda.min}, \code{converge.tol} and \code{nvar.max} are the options for stopping the fitting process. Smaller \code{lambda.min}, or smaller \code{converge.tol}, or larger \code{nvar.max} yields more accurate results, paricularly for deterministic computer experiments. \code{pen.norm} specifies the type of penalty norm in the loss function. \code{model} specifies the response type, which can be non-continuous response, in the case the loss function is replaced by negative log-likelihood function. More details can be seen in Sung et al. (2017+).
#'
#' @return An MRFA object is returned, for which \code{aic.MRFA}, \code{bic.MRFA} and \code{predict} methods exist.
#'
#' @seealso \code{\link{predict.MRFA}} for prediction of the MRFA model.
#' @author Chih-Li Sung <iamdfchile@gmail.com>
#'
#' @importFrom fields Wendland
#' @importFrom randtoolbox sobol
#' @importFrom plyr ldply
#' @importFrom utils combn
#' @import stats
#' @import grplasso
#' @import methods
#' @import foreach
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
#' @export
#'

MRFA_fit <- function(X, Y, weights = rep(1, length(Y)), order = 10, level = 10,
                     lambda.min = 1e-5, converge.tol = 1e-10, nvar.max = min(3 * length(Y), 3000), pen.norm = c("2", "N")[1],
                     model = LinReg(),
                     parallel = FALSE, verbose = TRUE){

  #########     setting      ########
  if(is.null(ncol(X))) X <- matrix(X, ncol = 1)
  n <- length(Y)
  if(nrow(X) != n) stop("Number of rows of X does not match the length of Y.")
  d <- ncol(X)
  min.x <- apply(X, 2, function(x) min(x))
  scale.x <- apply(X, 2, function(x) diff(range(x)))
  X <- t((t(X) - min.x)/scale.x)    # scale X to [0,1]
  if(order > d) order <- d

  #########     set gridpoints for basis functions     ########
  #########     set bandwidth for basis functions      ########
  Grid.ls <- generateGrid(order = order, level = level)
  gridpoint.ls <- Grid.ls$gridpoint.ls
  bandwidth.ls <- Grid.ls$bandwidth.ls

  #########     if native norm, precompute the cholesky decomposition     ########
  choldecompose.ls <- NULL
  if(pen.norm == "N"){
    choldecompose.ls <- vector("list", order)
    for(u in 1:order) choldecompose.ls[[u]] <- vector("list", level)

    for(u in 1:order){
      for(l in 1:level){
        grid.mx <- as.matrix(gridpoint.ls[[u]][[l]])
        group.ls <- list(list(effect = 1:u, resolution = l))
        Phi <- basis_fun(list(group.ls), grid.mx, gridpoint.ls, bandwidth.ls, parallel = parallel)
        choldecompose.ls[[u]][[l]] <- chol(Phi)
      }
    }
  }

  #########     set initial candidate group     ########
  candidate.group <- vector("list", d)
  for(j in 1:d) candidate.group[[j]][[1]] <- list(effect = j, resolution = 1)

  #########     set initial basis functions     ########
  Phi <- basis_fun(candidate.group, X, gridpoint.ls, bandwidth.ls, parallel)
  Phi <- cbind(rep(1,nrow(Phi)), Phi)

  group.index <- NA
  for(i in 1:length(candidate.group)){
    u <- length(candidate.group[[i]][[1]]$effect)
    l <- candidate.group[[i]][[1]]$resolution
    group.index <- c(group.index, rep(i, nrow(gridpoint.ls[[u]][[l]])))
  }

  #########     set lambda.max which allows only one group in the model    ########
  if (verbose)  cat("Setting lambda.max \n")
  lambda.test <- lambdamax(x = Phi, y = Y, index = group.index) + 100
  fit.tmp <- grplasso(x = Phi, y = Y, index = group.index, lambda = lambda.test, model = model,
                      control = grpl.control(trace = as.numeric(verbose)))
  sig.group <- unique(na.omit(group.index[fit.tmp$coef != 0]))
  while(length(sig.group) != 1){
    if(length(sig.group) > 1){
      lambda.test <- lambda.test + lambda.test/100
    }else{
      lambda.test <- lambda.test - lambda.test/100
    }
    fit.tmp <- grplasso(x = Phi, y = Y, index = group.index, lambda = lambda.test, model = model,
                        control = grpl.control(trace = as.numeric(verbose)))
    sig.group <- unique(na.omit(group.index[fit.tmp$coef != 0]))
  }
  lambda.max <- lambda.test
  coef.init <- c(fit.tmp$coef)
  if (verbose)  cat("Setting lambda.max =", lambda.max, " and lambda.min =", lambda.min, "\n")

  MRFA.fit <- grplasso_forward(x = Phi, y = Y, X = X, index = group.index, weights = weights,
                               order = order, level = level,
                               lambda.max = lambda.max, lambda.min = lambda.min, converge.tol = converge.tol, nvar.max = nvar.max,
                               pen.norm = pen.norm,
                               coef.init = coef.init, model = model,
                               candidate.group = candidate.group,
                               gridpoint.ls = gridpoint.ls,
                               bandwidth.ls = bandwidth.ls,
                               choldecompose.ls = choldecompose.ls,
                               penscale = sqrt, center = TRUE, standardize = TRUE,
                               control = grpl.control(trace = as.numeric(verbose)), parallel = parallel)

  MRFA.fit$min.x <- min.x
  MRFA.fit$scale.x <- scale.x
  invisible(MRFA.fit)
}
