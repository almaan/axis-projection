axis.projection <- function(object,
                            ab.column = "AB",
                            a.label = "A",
                            b.label ="B",
                            split.on = NULL){
  
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
  }
  
  
  idxs <- order(as.integer(unlist(idxs)))

  all.ab.projs <- do.call(rbind,all.ab.projs)
  all.ab.projs <- all.ab.projs[idxs,]
  
  return(Seurat::AddMetaData(object,all.ab.projs,col.name = new.col))
  
}


