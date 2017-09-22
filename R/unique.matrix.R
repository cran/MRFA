unique.matrix <- function(active.group, beta_hat, gridpoint.ls){
  active.group.unlist <- unlist(active.group, recursive = FALSE)
  ncol.vt <- sapply(active.group.unlist, function(x) nrow(gridpoint.ls[[length(x$effect)]][[x$resolution]]))
  cumsum.ncol <- cumsum(ncol.vt)

  unique.active.group <- lapply(active.group, function(x) x[length(x)])
  unique.ncol.vt <- sapply(unique.active.group, function(x) nrow(gridpoint.ls[[length(x[[1]]$effect)]][[x[[1]]$resolution]]))
  unique.beta <- rep(0, sum(unique.ncol.vt))
  cumsum.unique.ncol <- cumsum(unique.ncol.vt)

  for(ii in 1:length(unique.active.group)){
    group <- unique.active.group[[ii]][[1]]
    repeat.fg <- sapply(active.group.unlist, function(x) {
      if(length(x$effect) != length(group$effect) | x$resolution != group$resolution){
        return(FALSE)
      }else{
        all(x$effect == group$effect)
      }
    })

    plug.index <- cumsum.ncol[repeat.fg]
    beta_hat_sum <- rep(0, nrow(gridpoint.ls[[length(group$effect)]][[group$resolution]]))
    for(jj in 1:length(plug.index)){
      tmp.index <- (plug.index[jj] - nrow(gridpoint.ls[[length(group$effect)]][[group$resolution]]) + 1) : plug.index[jj]
      beta_hat_sum <- beta_hat_sum + beta_hat[tmp.index + 1] # plus an intercept
    }
    unique.beta[(cumsum.unique.ncol[ii] - length(beta_hat_sum) + 1):cumsum.unique.ncol[ii]] <- beta_hat_sum
  }

  unique.beta <- c(beta_hat[1], unique.beta)
  return(list(unique.beta = unique.beta, unique.active.group = unique.active.group))
}

