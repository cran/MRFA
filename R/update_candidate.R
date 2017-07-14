update_candidate <- function(add_active_group.index, active_group.index, candidate.group, order, level){
  ##### update candidate group according to strong effect heredity principle

  new_candidate.group <- list()
  kk <- 1
  for(group.index in add_active_group.index){
    active.group <- candidate.group[active_group.index]
    new_active.group <- candidate.group[[group.index]]

    ## update candidate group with increasing resolution
    new_active.resolution <- max(sapply(new_active.group, function(x) x$resolution))
    if(new_active.resolution < level){
      group.add <- new_active.group[sapply(new_active.group, FUN = function(x) x$resolution == new_active.resolution)]
      group.add <- lapply(group.add, function(x) {
        x$resolution <- x$resolution + 1
        return(x)})
      new_candidate.add <- c(new_active.group, group.add)

      new_candidate.group <- c(new_candidate.group, list(new_candidate.add))
    }

    ## update candidate group with interaction effect
    new_active.order <- max(sapply(new_active.group, function(x) length(x$effect)))
    if(new_active.order < order & length(active.group) > 0){
      same_order_level.group <- active.group[sapply(active.group, function(x.ls) max(sapply(x.ls, function(x) x$resolution)) == new_active.resolution &
                                                      max(sapply(x.ls, function(x) length(x$effect))) == new_active.order)]
      same_order_level.group <- c(same_order_level.group, list(new_active.group))
      if(length(same_order_level.group) > new_active.order){
        test.table <- combn(length(same_order_level.group) - 1, new_active.order)
        test.table <- rbind(test.table, length(same_order_level.group))
        group.tmp <- sapply(same_order_level.group, function(x.ls) x.ls[sapply(x.ls, function(x) length(x$effect) == new_active.order)])
        ## the effects which obey strong effect heredity principle
        effect_heredity.index <- which(apply(test.table, 2, function(x){
          group.tbl <- table(c(sapply(group.tmp[x], function(x.ls) x.ls$effect)))
          all(group.tbl == new_active.order) & length(group.tbl) == (new_active.order + 1)
        }))

        for(ii in effect_heredity.index){
          effect.new <- sort(unique(c(sapply(group.tmp[test.table[,ii]], function(x.ls) x.ls$effect))))
          group.new <- vector("list", new_active.resolution)
          for(jj in 1:new_active.resolution) group.new[[jj]] <- list(effect = effect.new, resolution = jj)
          new_candidate.add <- c(unique(unlist(same_order_level.group[test.table[,ii]], recursive = FALSE)), group.new)
          new_candidate.group <- c(new_candidate.group, list(new_candidate.add))
        }
      }
    }

    if(kk > 1) active_group.index <- c(active_group.index, group.index)  ## in case two active groups are entertained
  }

  if(length(new_candidate.group) > 1) new_candidate.group <- list.unique(new_candidate.group)
  if(length(new_candidate.group) > 0){
    delete.fg <- !is.anydupicate(new_candidate.group, candidate.group)
    if(any(delete.fg)){
      new_candidate.group <- new_candidate.group[delete.fg]
    }
  }
  return(new_candidate.group)
}
