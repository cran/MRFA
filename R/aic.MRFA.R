#' Extract AIC from a Fitted Multiresolution Functional ANOVA (MRFA) Model
#'
#' @description The function extracts Akaike information criterion (AIC) from a fitted MRFA model.
#'
#' @param fit a class MRFA object estimated by \code{MRFA_fit}.
#'
#' @return a vector with length \code{length(lambda)} returing AICs.
#'
#' @seealso \code{\link{predict.MRFA}} for prediction of the MRFA model.
#' @author Chih-Li Sung <iamdfchile@gmail.com>
#'
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
#' print(aic.MRFA(MRFA_model))
#' print(bic.MRFA(MRFA_model))
#'
#' #####   Prediction : AIC and BIC  ######
#' X.test <- matrix(runif(d*n_new), ncol = d)
#' Y.test <- apply(X.test, 1, otlcircuit)
#' lambda.aic <- MRFA_model$lambda[which.min(aic.MRFA(MRFA_model))]
#' Y.pred <- predict(MRFA_model, X.test, lambda = lambda.aic)$y_hat
#' print(sqrt(mean((Y.test - Y.pred)^2)))
#'
#' lambda.bic <- MRFA_model$lambda[which.min(bic.MRFA(MRFA_model))]
#' Y.pred <- predict(MRFA_model, X.test, lambda = lambda.bic)$y_hat
#' print(sqrt(mean((Y.test - Y.pred)^2)))
#' }
#' @export aic.MRFA
aic.MRFA <- function(fit)
{
  n <- length(fit$y)
  nloglik.vt <- fit$nloglik
  ## Update the negative log likelihood for linear regression
  if(fit$model@name == "Linear Regression Model"){
    nloglik.vt <- n * (log(2 * pi) + log(nloglik.vt) - log(n))/2
  }
  edf <- apply(fit$coefficients, 2, function(x) sum(x != 0))
  n <- length(fit$weights)
  2 * nloglik.vt + 2 * edf
}
