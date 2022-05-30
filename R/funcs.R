axis.projection <- function(object,
                            ab.column = "AB",
                            a.label = "A",
                            b.label ="B",
                            split.on = NULL,
                            include.orthogonal = FALSE){
  
  #' @export axis.projection
  
  new.col <- paste0(ab.column,"_projection")
  all.crd <- get.coordinates(object)
  
  
  idxs <- list()
  if (!(is.null(split.on))){
    split.labs <- unique(object@meta.data[[split.on]])
    for (sl in 1:length(split.labs)){
      idxs[[sl]] <- which(object@meta.data[[split.on]] == split.labs[sl])
    }
  } else {
    split.labs <- c(1)
    idxs[[1]] <- c(1:ncol(object))
  }

  
  all.ab.projs <- list()
  if (include.orthogonal){
    all.ab.orth.projs <- list()
  }
  
  
  for (sl in 1:length(split.labs)){
    crd <- all.crd[idxs[[sl]],]
    a.sel <- object@meta.data[idxs[[sl]],][[ab.column]] == a.label
    b.sel <- object@meta.data[idxs[[sl]],][[ab.column]] == b.label
    
    crd.a <-crd[a.sel,]
    crd.b <-crd[b.sel,]

    ab.vec <- get.axis.vector(crd.a,
                              crd.b)
    
    ab.proj <- project.on.vector(crd,
                                 ab.vec,
                                 reduce.dim = T,
                                 normalize = T)
    
    all.ab.projs[[sl]] <- ab.proj
    
    if (include.orthogonal){
      ab.orth.vec <- c(-ab.vec[2]/ab.vec[1],1)
      ab.orth.proj <-  project.on.vector(crd,
                                         ab.orth.vec,
                                         reduce.dim = T,
                                         normalize = T)
      all.ab.orth.projs[[sl]] <- ab.orth.proj
    }
  }
  
  
  idxs <- order(as.integer(unlist(idxs)))

  all.ab.projs <- do.call(rbind,all.ab.projs)
  all.ab.projs <- all.ab.projs[idxs,]
  new.meta.data <- data.frame(X = all.ab.projs)
  colnames(new.meta.data) <- new.col
  
  if (include.orthogonal){
    all.ab.orth.projs <- do.call(rbind,all.ab.orth.projs)
    all.ab.orth.projs <- all.ab.orth.projs[idxs,]
    new.meta.data[paste0("orth_to_",new.col)] <- all.ab.orth.projs
    
  }
  
  return(Seurat::AddMetaData(object,new.meta.data))
  
}

binner <- function(vals,
                   n.bins = 100,
                   lb.val = NULL,
                   ub.val = NULL
){
  
  #' @export binner
  
  
  if (is.null(lb.val)){
    mn <- min(vals)
  } else {
    mn <- lb.val
  }
  
  if (is.null(ub.val)){
   mx <- max(vals)
  } else {
    mx <- ub.val
  }
  
  delta <-  1/n.bins
  ss <- seq(mn,mx,delta)
  ns <- length(ss)
  binidx <-list()
  bin.av <- c()
  
  for (ii in 1:(ns-1)){
    up.val <- ss[ii+1]
    low.val <- ss[ii]
    in.bin <-which((vals >= low.val) & (vals < up.val))
    if (length(in.bin) > 0) {
      binidx[[toString(ii)]] <- in.bin
      bin.av <- c(bin.av,mean(vals[in.bin])) 
    }
  }
  
  names(binidx) <- c(1:length(binidx))
  bincenter <- ss + delta / 2
  bincenter <- bincenter[1:(length(bincenter)-1)]
  return(list(binval = bin.av,
              binidx = binidx
              )
        )
}

average.profiles <- function(X,
                             y,
                             n.bins = 100,
                             lb.val = NULL,
                             ub.val = NULL
                             ){
  #' @export average.profiles
  bin.data <- binner(y,n.bins = n.bins,lb.val = lb.val, ub.val = ub.val)
  new.y <- bin.data[["binval"]]
  new.n <- length(new.y) 
  new.X <- matrix(0,new.n,ncol = ncol(X))  
  colnames(new.X) <- colnames(X)
  rownames(new.X) <- paste0("bin_",c(1:new.n))

  for (ii in 1:new.n){
    idx <- bin.data[["binidx"]][[ii]]
    if (length(idx) > 1){
      new.X[ii,] <- colMeans(X[idx,]) 
    } else {
      new.X[ii,] <- X[idx,]
    }
  }
  return(list(X = new.X, y = new.y))
  
}
 
get.axis.genes <- function(X,
                           y,
                           bin.data = T,
                           n.bins = 100,
                           n.genes = NULL,
                           lb.val = NULL,
                           ub.val = NULL
                           ){
  
  #' @export get.axis.genes
  if (bin.data){
    binned.data <- average.profiles(X,
                                y,
                                n.bins = n.bins,
                                lb.val = lb.val, 
                                ub.val = ub.val)
    X.var <- binned.data$X
    y.var <- binned.data$y
  } else {
    X.var <- X
    y.var <- y
  }
  pen.fit <- penalized::penalized(y.var,
                       X.var,
                       lambda1=0,
                       lambda2 =1,
                       model="logistic")
  
  pen.order <- order(abs(pen.fit@penalized),
                     decreasing = T)
  
  
  n.genes <- ifelse(is.null(n.genes),length(pen.order),n.genes)
  axis.genes <- names(pen.fit@penalized)[pen.order[1:n.genes]]
  
  return(axis.genes)
  
}


