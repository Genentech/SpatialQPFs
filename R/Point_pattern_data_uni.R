#' A Point_pattern_data_uni function
#'
#' This is a function that calculates the features for univariate marker for spatial point pattern data. The methods included are Ripley's K-function, G-function and Clark and Evans nearest neighbor index
#'
#' @param path The path for the directory that contains the data
#' @param file The file name in the directory
#' @param cell_class_var The column name in each file that indicates the cell class variable
#' @param x_var The column name in each file that indicates the cell location x variable
#' @param y_var The column name in each file that indicates the cell location y variable
#' @param cell_type The cell type one wants to use as the univariate marker
#' @param scale The spatial range that user wants to investigate
#' @param myplot Whether to plot the results, if available, by default it is set as FALSE
#'
#' @return This function returns the univariate marker features for spatial point pattern data
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





Point_pattern_data_uni <- function(path = "/Users/lix233/Haystack/5862_cell_centers/",
                                        file = "f8008cf4-3fe7-4200-a2af-21c27708ae28.csv",
                                        cell_class_var = "cell_class",
                                        x_var = "X0",
                                        y_var = "X1",
                                        cell_type = "Lymphocyte", scale,
                                        myplot = FALSE){
  
  spp_df = read.csv(paste0(path, file)) %>% filter(get(cell_class_var) == cell_type)

  spp_df = spp_df %>% select("cell_class" = all_of(cell_class_var),
                             "x_coord" = all_of(x_var),
                             "y_coord" = all_of(y_var))

  L = max(max(spp_df$x_coord), max(spp_df$y_coord))

  spp_df = spp_df %>% mutate("x" = normalized_coords(spp_df$x_coord, L),
                             "y" = normalized_coords(spp_df$y_coord, L))


  radius = scale/L

  ##########################################################################################################################################
  ##########################################################################################################################################
  ####################################################        Spatial plots & features        ##############################################
  ##########################################################################################################################################
  ##########################################################################################################################################
  xy <- as.matrix(spp_df[, c('x','y')])

  bnd = owin(xrange = range(xy[,1]), yrange = range(xy[,2]))

  ln = with(spp_df,
            ppp(x = x, y = y, mark = cell_class, window = bnd) ## xrange and yrange
  )


  r = seq(0, radius, length.out = 20)

  ## G-function
  Gln = Gest(ln, r = r, correction = "km")
  if (myplot) {
    plot(Gln, lty = 1, main = paste0("G-function for ", cell_type))
  }  
  tmp = as.data.frame(cbind(Gln$km, Gln$r, Gln$theo)); colnames(tmp) = c('km', 'r', 'theo')
  tmp = tmp[is.finite(tmp$km), ]
  g_AUC = trapz(tmp$r, tmp$km-tmp$theo)
  g_r = tmp$r[which.max(tmp$km)]

  ## K-function
  Kln = Kest(ln, r = r, correction = "translate")
  if (myplot) {
    plot(Kln, lty = 1, main = paste0("Ripley's K-function for ", cell_type))
  }  
  tmp = as.data.frame(cbind(Kln$trans, Kln$r, Kln$theo)); colnames(tmp) = c('trans', 'r', 'theo')
  tmp = tmp[is.finite(tmp$trans), ]
  k_AUC = trapz(tmp$r, tmp$trans-tmp$theo)
  k_vals_med = quantile(tmp$trans, 0.5)
  k_vals_q1 = quantile(tmp$trans, 0.25)
  k_vals_q3 = quantile(tmp$trans, 0.75)
  k_vals_max = max(tmp$trans)

  ## Clark and Evans Aggregation Index
  CE = clarkevans(ln)[1]



  return(list("g_AUC" = g_AUC,
              "g_r" = g_r,
              "k_AUC" = k_AUC,
              "k_vals_med" = k_vals_med,
              "k_vals_q1" = k_vals_q1,
              "k_vals_q3" = k_vals_q3,
              "k_vals_max" = k_vals_max,
              "CE" = CE))

}
