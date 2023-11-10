post_sigma <- function (X, y, beta_hat, conf.level = 0.95,
                   var.estimation, w.estimation, K, nugget,
                   parallel, verbose)
{

  all.folds <- cv.folds(nrow(X), K)

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
    CIout <- foreach::foreach(k = seq(K), .combine = rbind, .verbose = verbose) %dopar%{
      omit <- all.folds[[k]]
      Phi.train <- X[-omit, , drop = FALSE]
      Phi.test <- X[omit, , drop = FALSE]
      Y <- y[-omit]
      Y.test <- y[omit]
      n <- nrow(Phi.train)
      n.xnew <- nrow(Phi.test)
      XtX <- t(Phi.train) %*% Phi.train
      Xsum <- matrix(apply(Phi.train, 2, sum), nrow = 1)
      res <- Y - Phi.train %*% t(beta_hat[,-1, drop=FALSE])
      U <- b <- rep(0, n.xnew)

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
        U[i] <- colSums((res - sum(a[-1] * beta_hat[,-1])) * c(res2  - sum(a[-1] * w)))
        b[i] <- sum(res2 - sum(a[-1] * w))
      }
      out <- cbind(U, b, Y.test)
    }
    U.vt <- CIout[,1]
    b.vt <- CIout[,2]
    Ytest.vt <- CIout[,3]
  }else{
    U.vt <- b.vt <- Ytest.vt <- rep(0, nrow(X))
    for (k in seq(K)) {
      omit <- all.folds[[k]]
      Phi.train <- X[-omit, , drop = FALSE]
      Phi.test <- X[omit, , drop = FALSE]
      Y <- y[-omit]
      Y.test <- y[omit]
      n <- nrow(Phi.train)
      n.xnew <- nrow(Phi.test)
      XtX <- t(Phi.train) %*% Phi.train
      Xsum <- matrix(apply(Phi.train, 2, sum), nrow = 1)
      res <- Y - Phi.train %*% t(beta_hat[,-1, drop=FALSE])
      U <- b <- rep(0, n.xnew)

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
        U[i] <- colSums((res - sum(a[-1] * beta_hat[,-1])) * c(res2  - sum(a[-1] * w)))
        b[i] <- sum(res2 - sum(a[-1] * w))
      }
      U.vt[omit] <- U
      b.vt[omit] <- b
      Ytest.vt[omit] <- Y.test
    }
  }

  froot <- function(x){
    c_alpha_1 <- U.vt + x * qnorm(1-(1-conf.level)/2) * sqrt(b.vt)
    c_alpha_2 <- U.vt - x * qnorm(1-(1-conf.level)/2) * sqrt(b.vt)

    UCI <- c_alpha_1/b.vt
    LCI <- c_alpha_2/b.vt
    sum(UCI > Ytest.vt & LCI < Ytest.vt) - conf.level * length(Ytest.vt)
  }

  sigma <- uniroot(froot, c(0, 20))$root

  return(sigma)

}
