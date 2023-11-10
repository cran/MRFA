grplasso_forward <- function(x, y, X, index, weights = rep(1, length(y)),
                             order = 10, level = 10,
                             lambda.max, lambda.min, converge.tol, nvar.max,
                             pen.norm = c("2", "N"),
                             coef.init = coef.init, model = LinReg(),
                             candidate.group = candidate.group,
                             gridpoint.ls = gridpoint.ls,
                             bandwidth.ls = bandwidth.ls,
                             choldecompose.ls = NULL, k = k,
                             penscale = sqrt, center = TRUE, standardize = TRUE, control = grpl.control(), parallel = FALSE,...)
{
  ###### modified from grplasso (see the original package grplasso by Lukas Meier)

  ## Set lambda based on lambda.max and lambda.min
  lambda_bandwidth <- 0.01
  if(lambda.min < 1e-5) lambda.min <- 1e-5
  lambda <- exp(seq(log(lambda.max), log(lambda.min), by = -lambda_bandwidth))

  ## Do some error checking
  ## Check the design matrix
  if(!is.matrix(x))
    stop("x has to be a matrix")

  if(any(is.na(x)))
    stop("Missing values in x not allowed!")

  ## Check the response
  if(!is.numeric(y))
    stop("y has to be of type 'numeric'")

  if(NROW(x) != length(y))
    stop("x and y have not correct dimensions")

  if(!model@check(y))
    stop("y has wrong format")

  ## Check the other arguments
  if(length(weights) != length(y))
    stop("length(weights) not equal length(y)")

  if(any(weights < 0))
    stop("Negative weights not allowed")

  if((!isTRUE(all.equal(weights, rep(weights[1], length(y))))) &
     (center | standardize))
    warning("Weights not considered for centering/scaling at the moment...")

  if(length(coef.init) != ncol(x))
    stop("length(coef.init) not equal ncol(x)")

  if(!is.numeric(index))
    stop("index has to be of type 'numeric'!")

  if(length(index) != ncol(x))
    stop("length(index) not equal ncol(x)!")

  if(is.unsorted(rev(lambda)))
    warning("lambda values should be sorted in decreasing order: lambda.min is not small enough.")

  if(all(is.na(index)))
    stop("None of the predictors are penalized.")

  check <- validObject(control) ## will stop the program if error occurs

  ## Extract the control information
  update.hess  <- control@update.hess
  update.every <- control@update.every
  inner.loops  <- control@inner.loops
  line.search  <- control@line.search
  max.iter     <- control@max.iter
  lower        <- control@lower
  upper        <- control@upper
  save.x       <- control@save.x
  save.y       <- control@save.y
  tol          <- control@tol
  trace        <- control@trace
  beta         <- control@beta
  sigma        <- control@sigma

  nrlambda <- length(lambda)
  ncolx    <- ncol(x)
  nrowx    <- nrow(x)

  if(nrlambda > 1 & update.hess == "always"){
    warning("More than one lambda value and update.hess = \"always\". You may want to use update.hess = \"lambda\"")
  }

  ## For the linear model, the Hessian is constant and has hence to be
  ## computed only *once*

  if(model@name == "Linear Regression Model"){
    if(update.hess != "lambda"){
      update.hess <- "lambda"
      if(trace >= 1)
        cat("Setting update.hess = 'lambda'\n")
    }
    if(update.every <= length(lambda)){
      update.every <- length(lambda) + 1
      if(trace >= 1)
        cat("Setting update.every = length(lambda) + 1\n")
    }
  }

  ## Which are the non-penalized parameters?
  any.notpen    <- any(is.na(index))
  inotpen.which <- which(is.na(index))
  nrnotpen      <- length(inotpen.which)

  intercept.which <- which(apply(x == 1, 2, all))
  has.intercept   <- length(intercept.which)

  if(!has.intercept & center){
    message("Couldn't find intercept. Setting center = FALSE.")
    center <- FALSE
  }

  if(length(intercept.which) > 1)
    stop("Multiple intercepts!")

  if(has.intercept)
    has.intercept.notpen <- is.na(index[intercept.which])
  else
    has.intercept.notpen <- FALSE

  others.notpen   <- nrnotpen - has.intercept.notpen
  notpen.int.only <- has.intercept.notpen & !others.notpen

  if(has.intercept & !center & standardize)
    warning("Are you sure that you don't want to perform centering in a model with intercept and standardized predictors?")

  ##if(center & others.notpen)
  if(others.notpen)
    warning("Penalization not adjusted to non-penalized predictors.")

  ## Index vector of the penalized parameter groups
  if(any.notpen){
    ipen <- index[-inotpen.which]
    ipen.which <- split((1:ncolx)[-inotpen.which], ipen)
  }else{
    if(has.intercept)
      warning("All groups are penalized, including the intercept.")
    ipen <- index
    ipen.which <- split((1:ncolx), ipen)
  }

  nrpen    <- length(ipen.which)
  dict.pen <- sort(unique(ipen))

  ## Table of degrees of freedom
  ipen.tab   <- table(ipen)[as.character(dict.pen)]

  x.old <- x

  if(center){
    if(!has.intercept) ## could be removed; already handled above
      stop("Need intercept term when using center = TRUE")

    mu.x                 <- apply(x[,-intercept.which], 2, mean)
    x[,-intercept.which] <- sweep(x[,-intercept.which], 2, mu.x)
  }

  ## Standardize the design matrix -> blockwise orthonormalization
  if(standardize){
    ##warning("...Using standardized design matrix.\n")
    stand        <- blockstand(x, ipen.which, inotpen.which)
    x            <- stand$x
    scale.pen    <- stand$scale.pen
    scale.notpen <- stand$scale.notpen
  }
  ## From now on x is the *normalized* design matrix!

  ## Extract the columns into lists, works faster for large matrices
  if(any.notpen){
    x.notpen <- list(); length(x.notpen) <- nrnotpen
    for(i in 1:length(inotpen.which))
      x.notpen[[i]] <- x[,inotpen.which[[i]], drop = FALSE]
  }

  x.pen <- list(); length(x.pen) <- length(nrpen)
  for(i in 1:length(ipen.which))
    x.pen[[i]] <- x[,ipen.which[[i]], drop = FALSE]

  Xjj <- list()

  ## Extract the needed functions
  check     <- validObject(model)
  invlink   <- model@invlink
  nloglik   <- model@nloglik
  ngradient <- model@ngradient
  nhessian  <- model@nhessian

  coef      <- coef.init
  coef.pen  <- coef.init
  if(any.notpen)
    coef.pen  <- coef[-inotpen.which]

  #norms.pen    <- c(sqrt(rowsum(coef.pen^2, group = ipen)))
  norms.pen <- rep(0, length(ipen.which))
  for(j in 1:length(ipen.which)) {
    ind <- ipen.which[[j]]
    if(any.notpen)
      ind <- ind - 1

    norms.pen[j] <- penalty_norm_fun(coef.pen[ind], candidate.group[[j]], choldecompose.ls, pen.norm = pen.norm)
  }

#   norms.pen.m  <- matrix(0, nrow = nrpen, ncol = nrlambda,
#                          dimnames = list(NULL, lambda))
#   norms.npen.m <- matrix(0, nrow = nrnotpen, ncol = nrlambda,
#                          dimnames = list(NULL, lambda))
  nloglik.v <- fn.val.v <- numeric(nrlambda)
  coef.m    <- grad.m   <- matrix(0, nrow = ncolx, ncol = nrlambda,
                                  dimnames = list(colnames(x), lambda))
#   fitted    <- linear.predictors <- matrix(0, nrow = nrowx, ncol = nrlambda,
#                                            dimnames = list(rownames(x), lambda))

  converged <- rep(TRUE, nrlambda)

  ## *Initial* vector of linear predictors (eta) and transformed to the
  ## scale of the response (mu)
  eta <- c(x %*% coef)
  mu <- invlink(eta)

  ## Create vectors for the Hessian approximations
  if(any.notpen){
    nH.notpen <- numeric(nrnotpen)
  }
  nH.pen <- numeric(nrpen)

  ## setting the initial values
  active_group.index <- vector(length = 0)
  add_group.check <- FALSE
  stop_loop.check <- FALSE

  ## start the loop
  pos <- 0
  #for(pos in 1:nrlambda){
  while(pos < nrlambda){
    pos <- pos + 1
    if(stop_loop.check) break
    l <- lambda[pos]
    l_lower <- l
    l_upper <- lambda[pos-1]
    if(pos == 1) size.updae <- 100
    else size.updae <- l_upper - l_lower
    lambda_update.check <- FALSE ## If one more active groups are entertained, then update lambda
    stop_update.check <- FALSE   ## If one more active groups are entertained, then update lambda and don't stop until add one more group.

    while(!stop_update.check){
      stop_update.check <- TRUE

      if(trace >= 2)
        cat("\nLambda:", l, "\n")

      ## Initial (or updated) Hessian Matrix of the *negative* log-likelihood
      ## function (uses parameter estimates based on the last penalty parameter
      ## value)

      if(update.hess == "lambda" & pos %% update.every == 0 | pos == 1){
        ## Non-penalized groups
        if(any.notpen){
          Xj <- x.notpen[[1]]
          nH.notpen[1] <- min(max(nhessian(Xj, mu, weights, ...), lower), upper)
        }
        ## Penalized groups
        for(j in 1:nrpen){
          ind <- ipen.which[[j]]
          Xj  <- x.pen[[j]]

          ## Save Xjj = t(X_j) %*% X_j for computation efficiency. only works for linear model
          if(model@name == "Linear Regression Model"){
            Xjj[[j]] <- crossprod(Xj, weights * Xj)
          }

          diagH <- numeric(length(ind))
          for(i in 1:length(ind)){
            diagH[i] <- nhessian(Xj[, i, drop = FALSE], mu, weights, ...)
          }
          nH.pen[j] <- min(max(diagH, lower), upper)
        }
      }

      ## Start the optimization process
      fn.val <- nloglik(y, eta, weights, ...) +
        l * sum(penscale(ipen.tab) * norms.pen)

      ## These are needed to get into the while loop the first time
      do.all <- FALSE
      d.fn   <- d.par <- 1

      counter    <- 1 ## Count the sub-loops
      iter.count <- 0 ## Count the loops through *all* groups
      fn.val.oldlambda <- fn.val

      ## Stop the following while loop if the convergence criterion is fulfilled
      ## but only if we have gone through all the coordinates

      ##while(d.fn > tol | d.par > sqrt(tol) | !do.all){
      while(d.fn > tol | d.par > sqrt(tol) | !do.all){
        ## Escape loop if maximal iteration reached
        if(iter.count >= max.iter){
          converged[pos] <- FALSE
          warning(paste("Maximal number of iterations reached for lambda[", pos,
                        "]", sep = ""))
          break
        }

        ## Save the parameter vector and the function value of the previous step
        fn.val.old <- fn.val
        coef.old   <- coef

        ## Check whether we have some useful information from the previous step

        ## Go through all groups if counter == 0 or if we have exceeded the
        ## number of inner loops (inner.loops)
        if(counter == 0 | counter > inner.loops){
          do.all <- TRUE
          guessed.active <- 1:nrpen
          counter <- 1
          if(trace >= 2)
            cat("...Running through all groups\n")
        }else{## Go through the groups which were identified at the previous step
          guessed.active <- which(norms.pen != 0)
          if(length(guessed.active) == 0){
            guessed.active <- 1:nrpen
            do.all <- TRUE
            if(trace >= 2)
              cat("...Running through all groups\n")
          }else{
            do.all <- FALSE
            if(counter == 1 & trace >= 2)
              cat("...Starting inner loop\n")
            counter <- counter + 1
          }
        }

        if(do.all)
          iter.count <- iter.count + 1

        ## These are used for the line search, start at initial value 1
        ## They are currently here for security reasons
        start.notpen <- rep(1, nrnotpen)
        start.pen    <- rep(1, nrpen)

        if(any.notpen){
          ## Optimize the *non-penalized* parameters
          j <- 1
          ind <- inotpen.which[j]
          Xj  <- x.notpen[[j]]

          ## Gradient of the negative log-likelihood function
          ngrad <- c(ngradient(Xj, y, mu, weights, ...))
          nH <- nH.notpen[j]

          ## Calculate the search direction
          d <- -(1 / nH) * ngrad
          ## Set to 0 if the value is very small compared to the current
          ## coefficient estimate

          d <- zapsmall(c(coef[ind], d), digits = 16)[2]


          ## If d != 0, we have to do a line search
          if(d != 0){
            scale <- min(start.notpen[j] / beta, 1) ##1
            coef.test      <- coef
            coef.test[ind] <- coef[ind] + scale * d

            Xjd       <- Xj * d
            eta.test  <- eta + Xjd * scale

            if(line.search){
              qh    <- sum(ngrad * d)

              fn.val0     <- nloglik(y, eta, weights, ...)
              fn.val.test <- nloglik(y, eta.test, weights, ...)

              qh <- zapsmall(c(qh, fn.val0), digits = 16)[1]

              ## Armijo line search. Stop if scale gets too small (10^-30).
              while(fn.val.test > fn.val0 + sigma * scale * qh & scale > 10^-30){
                ##cat("Doing line search (nonpen)\n")
                scale          <- scale * beta
                coef.test[ind] <- coef[ind] + scale * d
                eta.test       <- eta + Xjd * scale
                fn.val.test    <- nloglik(y, eta.test, weights, ...)
              } ## end if(line.search)
              if(scale <= 10^-30){ ## Do nothing in that case
                #cat("Running into problems with scale\n")
                #cat("qh", qh, "d", d, "\n")
                ## coef.test <- coef
                ## eta.test  <- eta
                ## mu        <- mu
                start.notpen[j] <- 1
              }else{ ## Update the information
                coef <- coef.test
                Xj  <- x.notpen[[j]]
                eta  <- eta.test
                eta.change1 <- TRUE
                mu   <- invlink(eta)
                start.notpen[j] <- scale
              }

              ## Save the scaling factor for the next iteration (in order that
              ## we only have to do very few line searches)
              ## start.notpen[j] <- scale

              ## Update the remaining information
              ## coef <- coef.test
              ## eta  <- eta.test
              ## mu   <- invlink(eta)
            } ## end if(abs(d) > sqrt(.Machine$double.eps))
          } ## end for(j in 1:nrnotpen)
        } ## if(any.notpen)

        ## Optimize the *penalized* parameter groups
        for(j in guessed.active){
          ind  <- ipen.which[[j]]
          npar <- ipen.tab[j]

          coef.ind       <- coef[ind]
          cross.coef.ind <- crossprod(coef.ind)

          ## Design matrix of the current group
          Xj <- x.pen[[j]]

          ## Negative gradient of the current group
          ngrad <- c(ngradient(Xj, y, mu, weights, ...))

          ## Update the Hessian if necessary
          if(update.hess == "always"){
            diagH <- numeric(length(ind))
            for(i in 1:length(ind)){ ## for loop seems to be faster than sapply
              diagH[i] <- nhessian(Xj[,i,drop = FALSE], mu, weights, ...)
            }
            nH <- min(max(diagH, lower), upper)
          }else{
            nH <- nH.pen[j]
          }

          cond       <- -ngrad + nH * coef.ind
          cond       <- c(cond) # v0.2
          cond.norm2 <- crossprod(cond)

          ## Check the condition whether the minimum is at the non-differentiable
          ## position (-coef.ind) via the condition on the subgradient.
          border <- penscale(npar) * l
          if(cond.norm2 > border^2){
            d <- (1 / nH) *
              (-ngrad - l * penscale(npar) * (cond / c(sqrt(cond.norm2))))
            ##d <- zapsmall(c(coef.ind, d))[-(1:npar)]
          }else{
            d <- -coef.ind
          }

          ## If !all(d == 0), we have to do a line search
          if(!all(d == 0)){
            scale <- min(start.pen[j] / beta, 1)

            coef.test      <- coef
            coef.test[ind] <- coef.ind + scale * d
            Xjd            <- c(Xj %*% d)
            eta.test       <- eta + Xjd * scale

            if(line.search){
              qh <- sum(ngrad * d) +
                l * penscale(npar) * sqrt(crossprod(coef.ind + d)) -
                l * penscale(npar)* sqrt(cross.coef.ind)

              if(model@name == "Linear Regression Model"){
                fn.val.test <- -2 * scale * crossprod(y - eta, weights * Xj) %*% d + scale^2 * t(d) %*% Xjj[[j]] %*% d # fn.val.test - fn.val0 (only for linear regression)
                fn.val0 <- 0

                left <- fn.val.test +
                  l  * penscale(npar) * sqrt(crossprod(coef.test[ind]))

                right <- fn.val0 + l  * penscale(npar) * sqrt(cross.coef.ind) +
                  sigma * scale * qh
              }else{
                fn.val.test    <- nloglik(y, eta.test, weights, ...)
                fn.val0        <- nloglik(y, eta, weights, ...)

                left <- fn.val.test +
                  l  * penscale(npar) * sqrt(crossprod(coef.test[ind]))

                right <- fn.val0 + l  * penscale(npar) * sqrt(cross.coef.ind) +
                  sigma * scale * qh
              }


              while(left > right & scale > 10^-30){
                ##cat("Doing line search (pen)\n")
                scale          <- scale * beta
                coef.test[ind] <- coef.ind + scale * d
                eta.test       <- eta + Xjd * scale
                if(model@name == "Linear Regression Model"){
                  fn.val.test <- -2 * scale * crossprod(y - eta, weights * Xj) %*% d + scale^2 * t(d) %*% Xjj[[j]] %*% d # fn.val.test - fn.val0 (only for linear regression)
                }else{
                  fn.val.test    <- nloglik(y, eta.test, weights, ...)
                }
                left <- fn.val.test +
                  l  * penscale(npar) * sqrt(crossprod(coef.test[ind]))

                right <- fn.val0 + l * penscale(npar) * sqrt(cross.coef.ind) +
                  sigma * scale * qh
              } ## end while(left > right & qh != 0)
            } ## end if(line.search)
            ## If we escaped the while loop because 'scale' is too small
            ## (= we add nothing), we just stay at the current solution to
            ## prevent tiny values
            if(scale <= 10^-30){ ## Do *nothing* in that case
              ##coef.test <- coef
              ##eta.test  <- eta
              ##mu        <- mu
              start.pen[j] <- 1
            }else{
              coef <- coef.test
              eta  <- eta.test
              mu   <- invlink(eta)
              start.pen[j] <- scale
            }
          } ## end if(!all(d == 0))
          #norms.pen[j] <- sqrt(crossprod(coef[ind]))
          norms.pen[j] <- penalty_norm_fun(coef[ind], candidate.group[[j]], choldecompose.ls, pen.norm = pen.norm)
        } ## end for(j in guessed.active)

        fn.val <- nloglik(y, eta, weights, ...) +
          l * sum(penscale(ipen.tab) * norms.pen)

        ## Relative difference with respect to parameter vector
        ##d.par <- sqrt(crossprod(coef - coef.old)) / (1 + sqrt(crossprod(coef)))

        d.par <-  max(abs(coef - coef.old) / (1 + abs(coef)))
        ##d.par <-  max(abs(coef - coef.old) / (ifelse(abs(coef), abs(coef), 1)))

        ## Relative difference with respect to function value (penalized
        ## likelihood)
        d.fn <- abs(fn.val.old - fn.val) / (1 + abs(fn.val))

        ## Print out improvement if desired (trace >= 2)
        if(trace >= 2){
          cat("d.fn:", d.fn, " d.par:", d.par,
              " nr.var:", sum(coef != 0), "\n")
        }

        ## If we are working on a sub-set of predictors and have converged
        ## we stop the optimization and will do a loop through all
        ## predictors in the next run. Therefore we set counter = 0.

        ##if(d.fn <= tol & d.par <= sqrt(tol)){
        if(d.fn <= tol & d.par <= sqrt(tol)){
          counter <- 0 ## will force a run through all groups
          if(trace >= 2 & !do.all)
            cat("...Subproblem (active set) solved\n")
        }
      } ## end of while(d.fn > tol | d.par > sqrt(tol) | !do.all)

      coef.m[,pos]            <- coef
#       fn.val.v[pos]           <- fn.val
#       norms.pen.m[,pos]       <- norms.pen
      nloglik.v[pos]          <- nloglik(y, eta, weights, ...)
#       linear.predictors[,pos] <- eta
#       fitted[,pos]            <- invlink(eta)

      if(trace == 1)
        cat("Lambda:", l, " nr.var:", sum(coef != 0), " ")

      ## update active groups
      all_active_group.index <- unique(na.omit(index[coef != 0]))
      # fix singularity on 12/24/18
      all_active_group.index <- sort(unique(c(all_active_group.index, active_group.index)))
      add_active_group.index <- all_active_group.index[!is.element(all_active_group.index, active_group.index)]

      ## If p > n then stop
      if(FALSE){ # close temperarily testing
        if(sum(coef != 0) > nrowx){
          active_group.index <- all_active_group.index
          coef.m            <- coef.m[,1:pos]
          #         fn.val.v          <- fn.val.v[1:pos]
          #         norms.pen.m       <- norms.pen.m[,1:pos]
          nloglik.v         <- nloglik.v[1:pos]
          #         linear.predictors <- linear.predictors[,1:pos]
          #         fitted            <- fitted[,1:pos]
          lambda            <- lambda[1:pos]
          stop_loop.check <- TRUE
          break
        }
      }

      ## Stop when it converges or number of active variables is greater than the pre-assigned number

      fn.diff <- abs(fn.val.oldlambda - fn.val) / (1 + abs(fn.val))
      if(trace == 1)
        cat("relative difference :", fn.diff, "\n")
      if(fn.diff < converge.tol | sum(coef != 0) > nvar.max){
        active_group.index <- all_active_group.index
        coef.m            <- coef.m[,1:pos, drop = FALSE]
        #         fn.val.v          <- fn.val.v[1:pos]
        #         norms.pen.m       <- norms.pen.m[,1:pos]
        nloglik.v         <- nloglik.v[1:pos]
        #         linear.predictors <- linear.predictors[,1:pos]
        #         fitted            <- fitted[,1:pos]
        lambda            <- lambda[1:pos]
        stop_loop.check <- TRUE

        if(trace == 1){
          cat("\n============= Entertaining new active groups: =============\n")
          print(lapply(candidate.group[add_active_group.index], function(x) x[[length(x)]]))
        }
        break
      }

      ## If lambda is increasing but no new active group in the current model (increase too much),
      ## then decrease lambda and update lambda

      if(lambda_update.check) {
        if(length(add_active_group.index) == 0){
          size.updae <- l - l_lower
          l_upper <- l
          l <- l - size.updae * 0.1
          if(size.updae/l_lower < 1e-7){
            l_upper <- l_lower
            l_lower <- max(l_lower - 1, (l_lower + lambda.min)/2)
            l <- l_lower
          }
          stop_update.check <- FALSE
          lambda <- c(lambda[1:(pos-1)], exp(seq(log(l), log(lambda.min), by = -lambda_bandwidth)))

          if(length(lambda) > nrlambda){  ## If the size of new lambda is greater than original lambda
            coef.m            <-   cbind(coef.m, matrix(0, ncol = length(lambda) - nrlambda, nrow = nrow(coef.m)))
#             fn.val.v          <-   c(fn.val.v, rep(0, length = length(lambda) - nrlambda))
#             norms.pen.m       <-   cbind(norms.pen.m, matrix(0, ncol = length(lambda) - nrlambda, nrow = nrow(norms.pen.m)))
            nloglik.v         <-   c(nloglik.v, rep(0, length = length(lambda) - nrlambda))
#             linear.predictors <-   cbind(linear.predictors, matrix(0, ncol = length(lambda) - nrlambda, nrow = nrow(linear.predictors)))
#             fitted            <-   cbind(fitted, matrix(0, ncol = length(lambda) - nrlambda, nrow = nrow(fitted)))
          }else if(length(lambda) < nrlambda){   ## If the size of new lambda is smaller than original lambda
            coef.m            <- coef.m[,1:length(lambda)]
#             fn.val.v          <- fn.val.v[1:length(lambda)]
#             norms.pen.m       <- norms.pen.m[,1:length(lambda)]
            nloglik.v         <- nloglik.v[1:length(lambda)]
#             linear.predictors <- linear.predictors[,1:length(lambda)]
#             fitted            <- fitted[,1:length(lambda)]
          }
          nrlambda <- length(lambda)
          next
        }
      }

      if(length(add_active_group.index) > 0 & pos < nrlambda){
        ## If more than one groups are entertained in current model, then increase lambda
        if(length(add_active_group.index) > 1){

          ## If lambda is not increased yet, then increase lambda to the previous lambda
          if(!lambda_update.check){
            l <- lambda[pos-1]
            lambda_update.check <- TRUE
            stop_update.check <- FALSE
            next
          }

          ## If lambda is already increased, increase lambda with a small increase step. Do this until one new active group is entertained.
          if(l != lambda[pos-1]){
            size.updae <- l_upper - l
            l_lower <- l
            l <- l + size.updae * 0.1
            if(size.updae/l_upper < 1e-7){
              l_lower <- l_upper
              l_upper <- l_upper + 1
              l <- l_upper
            }
            lambda <- c(lambda[1:(pos-1)], exp(seq(log(l), log(lambda.min), by = -lambda_bandwidth)))
            if(length(lambda) > nrlambda){
              coef.m            <-   cbind(coef.m, matrix(0, ncol = length(lambda) - nrlambda, nrow = nrow(coef.m)))
#               fn.val.v          <-   cbind(fn.val.v, matrix(0, ncol = length(lambda) - nrlambda, nrow = nrow(fn.val.v)))
#               norms.pen.m       <-   cbind(norms.pen.m, matrix(0, ncol = length(lambda) - nrlambda, nrow = nrow(norms.pen.m)))
              nloglik.v         <-   c(nloglik.v, rep(0, length = length(lambda) - nrlambda))
#               linear.predictors <-   cbind(linear.predictors, matrix(0, ncol = length(lambda) - nrlambda, nrow = nrow(linear.predictors)))
#               fitted            <-   cbind(fitted, matrix(0, ncol = length(lambda) - nrlambda, nrow = nrow(fitted)))
            }else if(length(lambda) < nrlambda){
              coef.m            <- coef.m[,1:length(lambda)]
#               fn.val.v          <- fn.val.v[,1:length(lambda)]
#               norms.pen.m       <- norms.pen.m[,1:length(lambda)]
              nloglik.v         <- nloglik.v[,1:length(lambda)]
#               linear.predictors <- linear.predictors[,1:length(lambda)]
#               fitted            <- fitted[,1:length(lambda)]
            }
            nrlambda <- length(lambda)
            lambda_update.check <- TRUE
            stop_update.check <- FALSE
            next
          }
        }

        ### Allow more than one new active groups in some cases
        #if(length(add_active_group.index) > 1) {
        #  warning("length(add_active_group.index) > 1 when lambda = ", l)
        #}

        if(trace == 1){
          cat("\n============= Entertaining new active groups: =============\n")
          print(lapply(candidate.group[add_active_group.index], function(x) x[[length(x)]]))
        }

        ## Update candidate group using strong effect heredity principle
        new_candidate.group <- update_candidate(add_active_group.index, active_group.index, candidate.group, order, level)
        active_group.index <- all_active_group.index
        if(length(new_candidate.group) == 0) next

        ## Update index based on the new candidate group
        index.old <- index
        index.max <- max(index, na.rm = TRUE)
        index.counter <- index.max + 1

        delete.index <- c()  ## in case p > n
        for(ii in 1:length(new_candidate.group)){
          index.add <- c()
          for(jj in 1:length(new_candidate.group[[ii]])){
            group.tmp <- new_candidate.group[[ii]][[jj]]
            index.add <- c(index.add, rep(index.counter, nrow(gridpoint.ls[[length(group.tmp$effect)]][[group.tmp$resolution]])))
          }
          ## If p > n then delete some of the candidate groups
          ## turn off for testing
          if(FALSE){
            if(length(index.add) <= nrow(x)){
              index <- c(index, index.add)
              index.counter <- index.counter + 1
            }else{
              delete.index = c(delete.index, ii)
            }
          }else{
            index <- c(index, index.add)
            index.counter <- index.counter + 1
          }
        }

        if(length(delete.index) > 0){
          new_candidate.group <- new_candidate.group[-delete.index]
        }

        if(length(new_candidate.group) == 0)  next

        if(trace == 1){
          cat("\n============= Adding new candiate groups: =============\n")
          print(lapply(new_candidate.group, function(x) x[[length(x)]]))
        }


        ## Update Phi
        candidate.group.old <- candidate.group
        candidate.group <- c(candidate.group, new_candidate.group)
        x.add <- try(basis_fun(new_candidate.group, X, gridpoint.ls, bandwidth.ls, k, parallel))

        x.pen.len <- length(x.pen)

        ## If memory is not enough to store new basis functions then stop and output results
        if(class(x.add)[1] == "try-error"){
          if(center){
            if(!has.intercept) ## could be removed; already handled above
              stop("Need intercept term when using center = TRUE")
            mu.x.old <- mu.x
          }
          index <- index.old
          candidate.group <- candidate.group.old
          coef.m            <- coef.m[,1:pos, drop=FALSE]
#           fn.val.v          <- fn.val.v[1:pos]
#           norms.pen.m       <- norms.pen.m[,1:pos]
          nloglik.v         <- nloglik.v[1:pos]
#           linear.predictors <- linear.predictors[,1:pos]
#           fitted            <- fitted[,1:pos]
          lambda            <- lambda[1:pos]
          stop_loop.check <- TRUE
          break
        }else{
          if(center){
            if(!has.intercept) ## could be removed; already handled above
              stop("Need intercept term when using center = TRUE")
            mu.x.old <- mu.x
            mu.x.add  <- apply(x.add, 2, mean)
            x.add     <- sweep(x.add, 2, mu.x.add)
            mu.x      <- c(mu.x, mu.x.add)
          }

          ## Standardize the design matrix -> blockwise orthonormalization
           ipen.which.add <- split(1:ncol(x.add), index[(ncolx + 1):length(index)])

          if(standardize){
            ##warning("...Using standardized design matrix.\n")
            scale.pen.old <- scale.pen
            stand        <- blockstand(x.add, ipen.which.add, vector(length = 0))
            x.add        <- stand$x
            scale.pen    <- c(scale.pen, stand$scale.pen)
          }

          ## From now on x is the *normalized* design matrix!
           for(i in 1:length(ipen.which.add))
             x.pen[[x.pen.len + i]] <- x.add[,ipen.which.add[[i]], drop = FALSE]
           col.add.n <- ncol(x.add)
        }

        ncolx <- ncolx + col.add.n
        if(any.notpen){
          ipen <- index[-inotpen.which]
          ipen.which <- split((1:ncolx)[-inotpen.which], ipen)
        }else{
          if(has.intercept)
            warning("All groups are penalized, including the intercept.")
          ipen <- index
          ipen.which <- split((1:ncolx), ipen)
        }

        nrpen.old <- nrpen
        nrpen    <- length(ipen.which)
        dict.pen <- sort(unique(ipen))

        ## Table of degrees of freedom
        ipen.tab   <- table(ipen)[as.character(dict.pen)]

        coef.m           <- rbind(coef.m, matrix(0, nrow = ncolx - nrow(coef.m), ncol = ncol(coef.m)))
        coef             <- c(coef, rep(0, ncolx - length(coef)))
        # norms.pen.m      <- rbind(norms.pen.m, matrix(0, nrow = length(ipen.which.add), ncol = ncol(norms.pen.m)))
        norms.pen        <- c(norms.pen, rep(0, length(ipen.which.add)))

        ## update
        for(j in (nrpen.old+1):nrpen){
          ind <- ipen.which[[j]]
          Xj  <- x.pen[[j]]

          if(model@name == "Linear Regression Model"){
            ## save Xjj = t(X_j) %*% X_j
            Xjj[[j]] <- crossprod(Xj, weights * Xj)
          }
          diagH <- numeric(length(ind))
          for(i in 1:length(ind)){
            diagH[i] <- nhessian(Xj[, i, drop = FALSE], mu, weights, ...)
          }
          nH.pen[j] <- min(max(diagH, lower), upper)
        }
        add_group.check <- TRUE
      }else{
        active_group.index <- all_active_group.index
        add_group.check <- FALSE
      }

    } ## end for(pos in 1:nrlambda){
  }

  ## Update column names since lambda is changed
  colnames(coef.m)  <- lambda
  #colnames(norms.pen.m) <-colnames(fitted) <- colnames(linear.predictors) <- lambda

  ## Transform the coefficients back to the original scale if the design
  ## matrix was standardized
  if(standardize){
    if(any.notpen)
      coef.m[inotpen.which,] <- (1 / scale.notpen) * coef.m[inotpen.which,]
    ## For df > 1 we have to use a matrix inversion to go back to the
    ## original scale
    for(j in 1:length(ipen.which)){
      ind <- ipen.which[[j]]
      coef.m[ind,] <- solve(scale.pen[[j]], coef.m[ind,,drop = FALSE])
    }
  }

  ## Need to adjust intercept if we have performed centering
  if(center){
    coef.m[intercept.which,] <- coef.m[intercept.which,] -
      apply(coef.m[-intercept.which,,drop = FALSE] * mu.x, 2, sum)
  }

  ## Overwrite values of x.old if we don't want to save it
  if(!save.x)
    x.old <- NULL
  if(!save.y)
    y <- NULL

  out <- list(x = x.old, ## use untransformed values
              y = y,
              order              = order,
              level              = level,
              lambda.min         = lambda.min,
              converge.tol       = converge.tol,
              converge.end       = fn.diff,
              nvar.max           = nvar.max,
              pen.norm           = pen.norm,
              gridpoint.ls       = gridpoint.ls,
              bandwidth.ls       = bandwidth.ls,
              coefficients       = coef.m,
              candidate.group    = candidate.group,
              active_group.index = active_group.index,
              active.group       = candidate.group[active_group.index],
              lambda             = lambda,
              index              = index,
              penscale           = penscale,
              model              = model,
              nloglik            = nloglik.v,
#               fitted             = fitted,
#               linear.predictors  = linear.predictors,
#               fn.val             = fn.val.v,
              converged          = converged,
              weights            = weights,
              control            = control,
              call               = match.call())
  structure(out, class = "MRFA")
}
