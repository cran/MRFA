is.anydupicate <- function(x1.ls, x2.ls){

  index.fg <- rep(FALSE, length(x1.ls))
  y1.ls <- lapply(x1.ls, function(x) x[[length(x)]])
  y2.ls <- lapply(x2.ls, function(x) x[[length(x)]])

  for(i in 1:length(x1.ls)){
    for(j in 1:length(x2.ls)){
      if(length(y1.ls[[i]]$effect) == length(y2.ls[[j]]$effect)){
        index.fg[i] <- all(y1.ls[[i]]$effect == y2.ls[[j]]$effect) &
          y1.ls[[i]]$resolution == y2.ls[[j]]$resolution
      }
    }
  }

  return(index.fg)
}
