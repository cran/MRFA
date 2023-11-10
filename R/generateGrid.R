generateGrid <- function(order, level){
#   #==================== old but good =========================================
#   gridpoint.ls[[3]]=list()
#   gridpoint.ls[[1]][[1]]=expand.grid(seq(0,1,0.2))
#   gridpoint.ls[[1]][[2]]=expand.grid(seq(0,1,0.1))
#   gridpoint.ls[[1]][[3]]=expand.grid(seq(0,1,0.05))
#
#   gridpoint.ls[[2]]=list()
#   gridpoint.ls[[2]][[1]]=expand.grid(seq(0,1,0.3),seq(0,1,0.3))
#   gridpoint.ls[[2]][[2]]=expand.grid(seq(0,1,0.2),seq(0,1,0.2))
#   gridpoint.ls[[2]][[3]]=expand.grid(seq(0,1,0.15),seq(0,1,0.15))
#
#   gridpoint.ls[[3]]=list()
#   gridpoint.ls[[3]][[1]]=expand.grid(seq(0,1,0.4),seq(0,1,0.4),seq(0,1,0.4))
#   gridpoint.ls[[3]][[2]]=expand.grid(seq(0,1,0.3),seq(0,1,0.3),seq(0,1,0.3))
#   gridpoint.ls[[3]][[3]]=expand.grid(seq(0,1,0.25),seq(0,1,0.25),seq(0,1,0.25))
#
#   dStar=0.75
#   bandwidth.ls=list()
#
#   #     bandwidths[[1]]=c(dStar,dStar/2,dStar/4)
#   #     bandwidths[[2]]=sqrt(2)*c(dStar,dStar/2,dStar/4)
#   #     bandwidths[[3]]=sqrt(3)*c(dStar,dStar/2,dStar/4)
#
#   bandwidth.ls[[1]]=c(dStar,dStar/2,dStar/4)
#   bandwidth.ls[[2]]=sqrt(2.5)*c(dStar,dStar/2,dStar/4)
#   bandwidth.ls[[3]]=sqrt(5)*c(dStar,dStar/2,dStar/4)

#   #==============================================================
#   gridpoint.ls <- vector("list", order)
#   for(u in 1:order) gridpoint.ls[[u]] <- vector("list", level)
#   bandwidth.ls <- vector("list", order)
#   for(u in 1:order) bandwidth.ls[[u]] <- rep(NA, level)
#
#   for(u in 1:min(order,4)){
#     for(l in 1:level){
#       if(l < 4){
#         ### gridpoint
#         grid <- seq(0, 1, length = max((6-u),2)*2^(l-1)+1)
#         grid.df <- data.frame(matrix(rep(grid, u), ncol = u))
#         gridpoint.ls[[u]][[l]] <- expand.grid(grid.df)
#
#         ### bandwidth
#         bandwidth.ls[[u]][l] <- 6 * 1/(length(grid)-1)
#       }else{
#         ### gridpoint
#         random.val <- sobol(n = nrow(gridpoint.ls[[u]][[l-1]]) * 2, dim = u)
#         gridpoint.ls[[u]][[l]] <- matrix(random.val, ncol = u)
#
#         ### bandwidth
#         n <- nrow(gridpoint.ls[[u]][[l]])
#         D <- gridpoint.ls[[u]][[l]]
#         Points <- expand.grid((1:n), (1:n))
#         Points <- cbind(Points[, 2], Points[, 1])
#         Points <- Points[Points[, 2] > Points[, 1], ]
#         junk <- (D[Points[, 2],,drop=FALSE] - D[Points[, 1],,drop=FALSE])^2
#         bandwidth.ls[[u]][l] <- max(3, (7-u)) * u *  min(sqrt(rowSums(junk))) / l
#       }
#     }
#   }
#
#   if(order > 4){
#     for(u in 5:order) {
#       for(l in 1:level){
#         ### gridpoint
#         random.val <- sobol(n = ((level+1)*level)/2  * max(3, (7-u))^u , dim = u)
#         random.val <- matrix(random.val, ncol = u)
#         gridpoint.ls[[u]][[l]] <- random.val[(((l-1)*l)/2*(7-u)^u+1):(((l+1)*l)/2*(7-u)^u),,drop=FALSE]
#
#         ### bandwidth
#         n <- length((((l-1)*l)/2*(7-u)^u+1):(((l+1)*l)/2*(7-u)^u))
#         D <- gridpoint.ls[[u]][[l]]
#         if(u == 1) D <- matrix(D, ncol = 1)
#         Points <- expand.grid((1:n), (1:n))
#         Points <- cbind(Points[, 2], Points[, 1])
#         Points <- Points[Points[, 2] > Points[, 1], ]
#         junk <- (D[Points[, 2],,drop=FALSE] - D[Points[, 1],,drop=FALSE])^2
#         bandwidth.ls[[u]][l] <- max(3, (7-u)) * u *  min(sqrt(rowSums(junk))) / l
#       }
#     }
#   }

  #====================== Final =================================
  #====================== sobol sequence ========================
  gridpoint.ls <- vector("list", order)
  bandwidth.ls <- vector("list", order)
  for(u in 1:order) {
    for(l in 1:level){
      random.val <- sobol(n = max(3, (7-u))^min(u,5) * l, dim = u)
      random.val <- matrix(random.val, ncol = u)
      ### gridpoint
      gridpoint.ls[[u]][[l]] <- random.val

      ### bandwidth
      if(l == 1){
        n <- nrow(random.val)
        D <- gridpoint.ls[[u]][[l]]
        if(u == 1) D <- matrix(D, ncol = 1)
        Points <- expand.grid((1:n), (1:n))
        Points <- cbind(Points[, 2], Points[, 1])
        Points <- Points[Points[, 2] > Points[, 1], ]
        junk <- (D[Points[, 2],,drop=FALSE] - D[Points[, 1],,drop=FALSE])^2
        bandwidth.ls[[u]][l] <-  6 * u * min(sqrt(rowSums(junk)))
      }else{
        bandwidth.ls[[u]][l] <- bandwidth.ls[[u]][1] / 2^(l-1)
      }
    }
  }

  for(u in 1:min(5, order)){
    if(u == 1){
      for(l in 1:level){
        ### gridpoint
        grid <- seq(0, 1, length = max((6-u),2)*2^(l-1)+1)
        grid.df <- data.frame(matrix(rep(grid, u), ncol = u))
        gridpoint.ls[[u]][[l]] <- expand.grid(grid.df)

        ### bandwidth
        bandwidth.ls[[u]][l] <- 6 * 1/(length(grid)-1)
      }
    }else{
      ### gridpoint
      grid <- seq(0, 1, length = max((6-u),2)+1)
      grid.df <- data.frame(matrix(rep(grid, u), ncol = u))
      gridpoint.ls[[u]][[1]] <- expand.grid(grid.df)

      ### bandwidth
      bandwidth.ls[[u]][1] <- 6 * 1/(length(grid)-1)
    }
  }


  #====================== good setting ========================
#   gridpoint.ls <- vector("list", order)
#   for(u in 1:order) gridpoint.ls[[u]] <- vector("list", level)
#   bandwidth.ls <- vector("list", order)
#   for(u in 1:order) bandwidth.ls[[u]] <- rep(NA, level)
#
#   for(l in 1:level) gridpoint.ls[[1]][[l]] <- expand.grid(seq(0, 1, length = 5*2^(l-1)+1))
#   for(l in 1:level) gridpoint.ls[[2]][[l]] <- expand.grid(seq(0, 1, length = l+4), seq(0, 1, length = l+4))
#
#   if(order > 2){
#     for(u in 3:order){
#       for(l in 1:level){
#         grid <- seq(0, 1, length = l+3)
#         grid.df <- data.frame(matrix(rep(grid, u), ncol = u))
#         gridpoint.ls[[u]][[l]] <- expand.grid(grid.df)
#       }
#     }
#   }
#
#   bandwidth.ls[[1]] <- 0.75 * 2^(1-1:level)
#   if(order > 1){
#     for(u in 2:order){
#       bandwidth.ls[[u]] <- 0.75 * sqrt(2.5 * (u-1)) * 2^(1-1:level)
#     }
#   }

#     #########     set gridpoints for basis functions     ########
#     gridpoint.ls <- vector("list", order)
#     for(u in 1:order) gridpoint.ls[[u]] <- vector("list", level)
#
#     for(u in 1:order){
#       for(l in 1:level){
#         spacing <- (u + 2^(2-l)-1)/10
#         grid <- seq(0, 1, spacing)
#         grid.df <- data.frame(matrix(rep(grid, u), ncol = u))
#         gridpoint.ls[[u]][[l]] <- expand.grid(grid.df)
#       }
#     }
#
#     #########     set bandwidth for basis functions     ########
#     bandwidth.ls <- vector("list", order)
#     for(u in 1:order) bandwidth.ls[[u]] <- rep(NA, level)
#     for(u in 1:order){
#       for(l in 1:level){
#         bandwidth.ls[[u]][l] <- sqrt((u^2+1)/2) * 0.75 / 2^(l-1)
#       }
#     }

  return(list(gridpoint.ls = gridpoint.ls, bandwidth.ls = bandwidth.ls))

}

