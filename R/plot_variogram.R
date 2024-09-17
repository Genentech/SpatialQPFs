#' Plot the semi-variogram/cross-variogram for geostatistics data
#'
#' This is a function that plot the semi-variogram/cross-variogram for geostatistics data
#'
#' @param v The sample variogram
#' @param m The fitted variogram 
#' @param title The title appears in the plot
#'
#' @return This function returns semi-variogram/cross-variogram plot for geostatistics data
#'
#' @import tidyverse
#' @import ggplot2
#' @import gstat
#'
#' @author Xiao Li, \email{xiao.li.xl2@roche.com}
#'
#' @export



plot_variogram <- function(v, m, title) {
  preds = variogramLine(m, maxdist = max(v$dist))
  print(ggplot() + 
          geom_point(data = v, aes(x = dist, y = gamma, size=np)) +
          geom_line(data = preds, aes(x = dist, y = gamma)) +
          ggtitle(title) +
          theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold")))
}