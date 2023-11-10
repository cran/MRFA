list.unique <- function(x.ls){
  index.fg <- rep(TRUE, length(x.ls))
  y.ls <- lapply(x.ls, function(x) x[[length(x)]])
  for(i in 1:(length(x.ls)-1)){
    for(j in (i+1):length(x.ls)){
      if(length(y.ls[[i]]$effect) == length(y.ls[[j]]$effect)){
        index.fg[j] <- !(all(y.ls[[i]]$effect == y.ls[[j]]$effect) &
          y.ls[[i]]$resolution == y.ls[[j]]$resolution)
      }
    }
  }
  x.ls[index.fg]
}
