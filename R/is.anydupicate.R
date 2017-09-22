is.anydupicate <- function(x1.ls, x2.ls, type = c(1,2)[1]){

  if(type == 1){
    index.fg <- rep(FALSE, length(x1.ls))
    y1.ls <- lapply(x1.ls, function(x) x[[length(x)]])
    y2.ls <- lapply(x2.ls, function(x) x[[length(x)]])

    for(i in 1:length(x1.ls)){
      index.fg.tmp <- rep(FALSE, length(x2.ls))
      for(j in 1:length(x2.ls)){
        if(length(y1.ls[[i]]$effect) == length(y2.ls[[j]]$effect)){
          index.fg.tmp[j] <- all(y1.ls[[i]]$effect == y2.ls[[j]]$effect) &
            y1.ls[[i]]$resolution == y2.ls[[j]]$resolution
        }
        index.fg[i] <- any(index.fg.tmp)
      }
    }
  }else{
    index.fg <- rep(FALSE, length(x1.ls))
    y1.ls <- x1.ls
    y2.ls <- x2.ls
    for(i in 1:length(x1.ls)){
      index.fg.tmp <- rep(FALSE, length(x2.ls))
      for(j in 1:length(x2.ls)){
        if(length(y1.ls[[i]]$effect) == length(y2.ls[[j]]$effect)){
          index.fg.tmp[j] <- all(y1.ls[[i]]$effect == y2.ls[[j]]$effect) &
            y1.ls[[i]]$resolution == y2.ls[[j]]$resolution
        }
      }
      index.fg[i] <- any(index.fg.tmp)
    }
  }

  return(index.fg)
}
