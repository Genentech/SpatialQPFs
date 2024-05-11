#' Data visualization
#'
#' This is a function that visualizes the input cell level data, including raw spatial map and smoothed density map
#'
#' @param path The path for the directory that contains the data
#' @param file The file name in the directory
#' @param cell_class_var The column name in each file that indicates the cell class variable
#' @param x_var The column name in each file that indicates the cell location x variable
#' @param y_var The column name in each file that indicates the cell location y variable
#' @param cell_type The cell type one wants to use as the univariate marker
#'
#' @return This function returns spatial distribution plot for certain cell type
#'
#' @import tidyverse
#' @import spatstat
#' @import plotly
#'
#' @author Xiao Li, \email{li.xiao@gene.com}
#'
#' @export





Data_Vis <- function(path = "/Users/lix233/Haystack/5862_cell_centers/",
                     file = "f8008cf4-3fe7-4200-a2af-21c27708ae28.csv",
                     cell_class_var = "cell_class",
                     x_var = "X0",
                     y_var = "X1",
                     cell_type = "Lymphocyte"){
  
  spp_df = read.csv(paste0(path, file)) %>% filter(get(cell_class_var) == cell_type)

  spp_df = spp_df %>% select("cell_class" = all_of(cell_class_var),
                             "x" = all_of(x_var),
                             "y" = all_of(y_var))

  
  p1 = ggplot(spp_df, aes(x = x, y = y, col = cell_class)) + geom_point(colour = "red") +
    xlab('') + ylab('') + guides(col = FALSE) + ggtitle(paste0("Spatial distribution of ", cell_type)) +
    theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
          strip.text = element_text(size = 15, face = "bold")) + 
    guides(x = "none", y = "none")
  print(p1)
  
  
  
  p2 = ggplot(spp_df, aes(x = x, y = y)) + 
    coord_equal() + 
    xlab('') + 
    ylab('') + 
    stat_density2d(aes(fill = ..level..),
                   geom = "polygon") + 
    scale_fill_viridis_c() + 
    theme(legend.position = 'none') + 
    guides(x = "none", y = "none") + 
    ggtitle(paste0("Spatial density of ", cell_type)) +
    theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"))
  print(p2)
  
  
  
  xy <- as.matrix(spp_df[, c('x','y')])
  
  bnd = owin(xrange = c(0, max(xy[,1])), yrange = c(0, max(xy[,2])))
  
  
  ln = with(spp_df,
            ppp(x = x, y = y, marks = cell_class, window = bnd)
  )
  
  d = data.frame(density(subset(ln, marks == cell_type), edge=TRUE, diggle=TRUE))
  p3 = plot_ly(d, x = ~x, y = ~y, z = ~value)
  
  print(p3)
  

}
