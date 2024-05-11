#' Discrimination of target cell population
#'
#' This is a function that discriminates a target cell population into 2 subgroups, in spatial relationships with the reference cell population. 
#' In case of lymphocytes and tumor cells, the lymphocytes are discriminated into intra-tumoral lymphocytes (ITL) and adjacent-tumoral lymphocytes (ATL)
#'
#' @param path The path for the directory that contains the data
#' @param file The file name in the directory
#' @param cell_class_var The column name in each file that indicates the cell class variable
#' @param x_var The column name in each file that indicates the cell location x variable
#' @param y_var The column name in each file that indicates the cell location y variable
#' @param from_type The cell type one wants to use as the "from" cell type
#' @param to_type The cell type one wants to use as the "to" cell type
#' @param micron_per_pixel The ratio of pixel to micron, i.e. scan resolution
#' @param myplot Whether to plot the results, if available, by default it is set as FALSE
#'
#' @return This function returns the ratio of ITL vs. target cell population, and the ratio of ITL vs. reference cell population
#'
#' @import tidyverse
#' @import pracma
#' @import splancs
#' @import rsdepth
#' @import mclust
#' @import FNN
#' @import spatstat
#' @import polyCub
#' @import dbmss
#' @import ecespa
#' @import spdep
#' @import gstat
#'
#' @author Xiao Li, \email{li.xiao@gene.com}
#'
#' @export







Point_pattern_data_ITLR <- function(path = "/Users/lix233/Haystack/5862_cell_centers/",
                                         file = "f8008cf4-3fe7-4200-a2af-21c27708ae28.csv",
                                         cell_class_var = "cell_class",
                                         x_var = "X0",
                                         y_var = "X1",
                                         from_type = "Lymphocyte",
                                         to_type = "Tumor",
                                         micron_per_pixel = 1,
                                         myplot = FALSE){

  spp_df = read.csv(paste0(path, file))

  spp_df = spp_df %>% select("cell_class" = all_of(cell_class_var),
                             "x_coord" = all_of(x_var),
                             "y_coord" = all_of(y_var))


  ## change from pixel to micron for computation mem issue
  spp_df$x_coord = round(spp_df$x_coord*micron_per_pixel)
  spp_df$y_coord = round(spp_df$y_coord*micron_per_pixel)

  spp_df = spp_df %>% select(row = y_coord, col = x_coord, cell_class)


  spp_df$row = spp_df$row - min(spp_df$row) + 1
  spp_df$col = spp_df$col - min(spp_df$col) + 1


  cell.c = spp_df %>% filter(cell_class == to_type) %>% select(x = col, y = row)


  cv.tumor = mse2d(as.points(cell.c), poly = cbind(c(min(cell.c$x), min(cell.c$x), max(cell.c$x), max(cell.c$x)),
                                                   c(min(cell.c$y), max(cell.c$y), max(cell.c$y), min(cell.c$y))),
                   nsmse = 10 , range = 0.1*min(max(cell.c$x)-min(cell.c$x), max(cell.c$y)-min(cell.c$y)))

  h = cv.tumor$h[which.min(cv.tumor$mse)]
    
  cell.l = spp_df %>% filter(cell_class == from_type) %>% select(x = col, y = row)

  nx = length(seq(min(spp_df$col), max(spp_df$col), 2))
  ny = length(seq(min(spp_df$row), max(spp_df$row), 2))

  tumor_density = kernel2d(as.points(cell.c), poly = cbind(c(min(cell.c$x), min(cell.c$x), max(cell.c$x), max(cell.c$x)),
                                                           c(min(cell.c$y), max(cell.c$y), max(cell.c$y), min(cell.c$y))),
                           h0 = h, nx = nx, ny = ny)
  spp_density = setNames(reshape2::melt(tumor_density$z), c('x_coords', 'y_coords', 'density'))
  x_cols = cbind(1:nx, seq(min(spp_df$col), max(spp_df$col), length.out = nx)); colnames(x_cols) = c('x_coords', 'x')
  y_rows = cbind(1:ny, seq(min(spp_df$row), max(spp_df$row), length.out = ny)); colnames(y_rows) = c('y_coords', 'y')

  left_join(spp_density, as.data.frame(x_cols), by = 'x_coords') -> temp
  temp <- left_join(temp, as.data.frame(y_rows), by = 'y_coords')
  sp_prox = rep(NA, nrow(cell.l))
  res = get.knnx(data = temp[,c('x','y')], query = cell.l[,c('x','y')], k=1)[[1]]
  for (i in 1:nrow(cell.l)) {
    id = res[i,]
    sp_prox[i] = temp[as.numeric(id), 'density']
  }



  if (length(unique(sp_prox)) > 2){
    fit = Mclust(sp_prox, G = 2)
    immune_class = fit$classification
    immune_class[immune_class == which.max(fit$parameters$mean)] = 'ITL'
    immune_class[immune_class == which.min(fit$parameters$mean)] = 'ATL'
    ITLR = sum(immune_class == 'ITL')/length(immune_class)
    ITLR2 = sum(immune_class == 'ITL')/nrow(cell.c)
  } else if (length(unique(sp_prox)) == 2) {
    immune_class = rep('ATL', length(sp_prox))
    immune_class[sp_prox == max(unique(sp_prox))] = "ITL"
    ITLR = sum(immune_class == 'ITL')/length(immune_class)
    ITLR2 = sum(immune_class == 'ITL')/nrow(cell.c)
  } else {
    ITLR = ITLR2 = 0
  }

  if(myplot){
    plot(cell.c[,1], cell.c[,2], col = "black", pch=1, cex=1, xlab='x', ylab='y', main = "ATL (green), ITL (red), Reference cells (black)")
    points(cell.l[immune_class=='ITL', 'x'], cell.l[immune_class=='ITL', 'y'], pch=1, cex=1, col='red')
    points(cell.l[immune_class=='ATL', 'x'], cell.l[immune_class=='ATL', 'y'], pch=1, cex=1, col='green')
  }
  

  return(list("ITLR" = ITLR, "ITLR2" = ITLR2))
}
