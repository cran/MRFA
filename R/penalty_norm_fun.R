penalty_norm_fun <- function(beta, guess.group, choldecompose.ls, pen.norm = c("2", "N")){

  if(pen.norm == "2" | pen.norm == 2){
    norms.pen <- sqrt(crossprod(beta))
  }else{
    count.index <- 0
    norms.pen <- 0
    for(i in 1:length(guess.group)){
      group <- guess.group[[i]]
      u <- group$effect
      l <- group$resolution
      L <- choldecompose.ls[[length(u)]][[l]]

      norms.pen <- norms.pen + crossprod(L %*% beta[(count.index + 1):(count.index + nrow(L))])
      if(is.na(norms.pen)) stop()
      count.index <- count.index + nrow(L)
    }
    norms.pen <- sqrt(norms.pen)
  }
  return(norms.pen)
}
