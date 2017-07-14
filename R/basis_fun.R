basis_fun <- function(group.ls, X, gridpoint.ls, bandwidth.ls, parallel = FALSE){
  ### Construct basis function based on candidate group, gridpoints and bandwidth

  group.ls <- unlist(group.ls, recursive = FALSE)

  Phi <- ldply(group.ls, function(xx){
    Dist <- t(apply(gridpoint.ls[[ length(xx$effect) ]][[ xx$resolution ]], 1, function(node){
      X.u <- matrix(X[,xx$effect], ncol = length(xx$effect))
      dist <- sweep(X.u, 2, node, FUN = "-")
      bandwith <- bandwidth.ls[[ length(xx$effect) ]][ xx$resolution ]
      apply(dist, 1, FUN = function(x) sqrt(sum(x^2))/bandwith)
    }))
    if(nrow(X) == 1) Dist <- t(Dist)   # for only one testing data case
    Wendland(d = Dist, dimension = length(xx$effect), k = 2)
  }, .parallel = parallel)

  return(t(Phi))
}
