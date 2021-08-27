plot.axis.projection <- function(object,
                                 features,
                                 ab.column = "AB",
                                 a.label = "A",
                                 b.label = "B",
                                 split.on = NULL,
                                 span = 0.3,
                                 normalize = F
){
  
  features <- c(features) 
  n.features <- length(features)
  
  plot.df <- as.data.frame(t(as.matrix(Seurat::GetAssay(object)[features,])))
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
  
  fill_color <- viridis::scale_fill_viridis(discrete = T,option="magma")
  color_color <- viridis::scale_color_viridis(discrete = T,option="viridis")
  
  g <- ggplot2::ggplot(data = plot.df,
              aes_string(x = "x",
                         y="value",
                         fill=ifelse(is.null(split.on),"feature",split.on),
                         color="feature"
                        )) +
    ggplot2::geom_point(size = 1, alpha = 0.1)+
    ggplot2::geom_smooth(method="loess",
                alpha=0.3,
                span = span) +
    ggplot2::xlim(c(0,1)) +
    ggplot2::scale_x_continuous(paste0(a.label,"-",b.label," axis"),
                       breaks = c(0,1),
                       labels = c(a.label,b.label)) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(),
          axis.line = ggplot2::element_line(colour = "black"),
          axis.text.x = ggplot2::element_text(size = 20),
          axis.text.y = ggplot2::element_text(size = 15),
          axis.title = ggplot2::element_text(size=14,
                                  face="bold")) +
    ggplot2::ylab("Feature Value") +
    ggplot2::ggtitle("Axis projection") +
    fill_color +
    color_color 
  
  return(g)
}