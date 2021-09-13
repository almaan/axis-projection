plot.axis.projection <- function(object,
                                 features,
                                 ab.column = "AB",
                                 a.label = "A",
                                 b.label = "B",
                                 split.on = NULL,
                                 span = 0.3,
                                 show.data = T,
                                 title = NULL,
                                 fill.color = NULL,
                                 color.color = NULL,
                                 data.alpha = 0.1,
                                 fill.alpha = 0.3,
                                 xlabel = NULL,
                                 ylabel = NULL,
                                 normalize = F
){
  
  #' @export plot.axis.projection
  
  features <- c(features) 
  n.features <- length(features)
  
  plot.df <- list()
  
  for (feature in features){
    if (feature %in% rownames(object)){
      tmp.feature <- as.numeric(Seurat::GetAssay(object)[feature,])
    } else if (feature %in% colnames(object@meta.data)){
      tmp.feature <- object@meta.data[feature]
    } else {
      print(sprintf("ERORR %s not in meta or assay data",feature))
    }
    plot.df[[feature]] <- tmp.feature
  }
    
  plot.df <- as.data.frame(do.call(cbind,plot.df))
  if (normalize){
    plot.df <- plot.df / apply(plot.df,2,max)
  }
  
  
  plot.df["x"] <- object@meta.data[paste0(ab.column,"_projection")]
  if (!(is.null(split.on))){
    plot.df[split.on] <- object@meta.data[split.on]
    id.vars <- c("x",split.on)
  } else {
    id.vars <- c("x")
  }
  plot.df <- reshape2::melt(plot.df,
                  id.vars = id.vars,
                  variable.name = "feature")
  
  if (is.null(fill.color)){
    fill.color <- viridis::scale_fill_viridis(discrete = T,option="magma")
  }
  
  if (is.null(color.color)){
    color.color <- viridis::scale_color_viridis(discrete = T,option="viridis")
  }

  g <- ggplot2::ggplot(data = plot.df,
              ggplot2::aes_string(x = "x",
                         y="value",
                         fill=ifelse(is.null(split.on),
                                     "feature",split.on),
                         color="feature"
                        ))
  if (show.data){
    g <- g + ggplot2::geom_point(size = 1, alpha = data.alpha)
  }
    
    g <- g + ggplot2::geom_smooth(method="loess",
                alpha=fill.alpha,
                span = span) +
    ggplot2::scale_x_continuous(ifnull(xlabel,paste0(a.label,"-",b.label," axis")),
                       breaks = c(0,1),
                       labels = c(a.label,b.label),
                       limits = c(0,1)
                       ) +
    
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(),
          axis.line = ggplot2::element_line(colour = "black"),
          axis.text.x = ggplot2::element_text(size = 20),
          axis.text.y = ggplot2::element_text(size = 15),
          axis.title = ggplot2::element_text(size=14,
                                  face="bold")) +
    
    ggplot2::ylab(ifnull(ylabel,"Feature Value")) +
    ggplot2::ggtitle(ifnull(title,"Axis projection")) +
    fill.color +
    color.color 
  
  return(g)
}

