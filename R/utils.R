get.coordinates <- function(object){
  return(object@tools$Staffli@meta.data[c("pixel_x","pixel_y")])
}

l2.norm <- function(x){
  return(sqrt(sum(x**2)))
}

l2.distance <- function(x,y){
  return(l2.norm(x-y))
}

project.on.vector <- function(data,v,reduce.dim = FALSE,normalize = TRUE){
  nv <- v / l2.norm(v)
  xv <- as.matrix(data) %*% nv
  if (reduce.dim){
    pr<- xv
    if (normalize){
      mn <- min(pr)
      mx <- max(pr)
      pr <- (pr - mn) / (mx-mn)
    }
  } else {
    pr <- t(as.matrix(drop(outer(nv,xv))))
    if (normalize) {
      cmn <- apply(pr,2,min)
      cmx <- apply(pr,2,max)
      pr <- (pr - cmn) / (cmx-cmn)
    }
  }
  return(pr)
}

get.axis.vector <- function(a.crd,
                            b.crd,
                            n.max = 200){
  
  na <- nrow(a.crd)
  nb <- nrow(b.crd)
  if (nrow(a.crd) > n.max){
    x.a.crd <- a.crd[sample(na,n.max,replace = F),]
  } else {
    x.a.crd <- a.crd
  }
  
  if (nrow(b.crd) > n.max){
    x.b.crd <- b.crd[sample(nb,n.max,replace = F),]
  } else {
    x.b.crd <- b.crd
  }
  
  pw.diff <- nice.diff.vector(x.a.crd,x.b.crd)
  pw.diff <- pw.diff / apply(pw.diff,1,l2.norm)
  axis.vector <- colMeans(pw.diff)
  return(axis.vector)
  
}

nice.diff.vector <- function(x,y){
  x <- as.matrix(x)
  y <- as.matrix(y)
  nx <- nrow(x)
  ny <- nrow(y)
  res <- matrix(0,nx,2)
  
  for (ii in c(1:nx)){
      idx.y <-which.min(apply(y,1,l2.distance,y=x[ii,]))
      res[ii,] <- y[idx.y,] -  x[ii,]
  }    
  return(res)
}


