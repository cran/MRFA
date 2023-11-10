plotCVMRFA <- function(object){
  lambda <- object$lambda
  cv <- object$cv
  cv.error <- object$cv.error
  plot(lambda, cv, type = "b", ylim = range(cv, cv + cv.error, cv - cv.error), xlab = "lambda", ylab="Cross-Validated MSE")
  error.bars(lambda, cv + cv.error, cv - cv.error,
             width = 1/length(lambda))
  invisible()
}


