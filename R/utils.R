get.coordinates <- function(object,xcrd = "pixel_x",ycrd = "pixel_y"){
  #' @export get.coordinates
  return(object@tools$Staffli@meta.data[c(xcrd,ycrd)])
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
  #' @export nice.diff.vector
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


ifnull <- function(x,y){
  #' @export ifnull
  return(ifelse(is.null(x),y,x))
}


trapz <- function(x,y,return.sum = T,interpolate = T,delta = 0.1){
  n <- length(x)
  if (interpolate){
    x.mod <- seq(x[1],x[n],delta)
    y.mod <- pracma::interp1(x,y,x.mod)
  } else {
    x.mod <- x
    y.mod <- y
    delta <- diff(x)
  }
  n <- length(x.mod)
  a.raw <- (y.mod[1:(n-1)] + y.mod[2:n]) * delta / 2
  a.cs <- cumsum(a.raw)
  if (return.sum){
    return(list(x = x.mod, y= y.mod,area = a.cs[n-1]))
  } else{
    return(list(x = x.mod, y = y.mod, area = a.cs))
  }
}

get.mass.cutoff <- function(x,y,p = 0.5,interpolate = T, delta = 0.1){
  int.res <- trapz(x,y,return.sum = F,interpolate = interpolate, delta = delta)
  a.raw <- int.res[["area"]]
  a.fraction <- a.raw / a.raw[length(a.raw)]
  a.diff <- abs(a.fraction - p)
  val <- int.res[["x"]][which(a.diff == min(a.diff))]
  if (length(val) > 1){
    val <- val[1]
  }
  return(val)
}
