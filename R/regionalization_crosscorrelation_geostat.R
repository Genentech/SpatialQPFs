#' Cross-variogram for geostatistics data
#'
#' This is the core function that calculates the features for cross-variogram, using geostatistics methods
#'
#' @param spp_df The input spatial data.frame, need to have 3 columns: "x_coord", "y_coord" and "cell_class"
#' @param from_type The cell type one wants to use as the "from" cell type
#' @param to_type The cell type one wants to use as the "to" cell type
#' @param scale The spatial range that user wants to investigate
#' @param myplot Whether to plot the results, if available, by default it is set as FALSE
#'
#' @return This core function calculates the features for cross-variogram
#'
#'
#' @import magrittr
#' @import dplyr
#' @import ggplot2
#' @importFrom gstat variogram fit.variogram variogramLine vgm
#'
#' @author Xiao Li, \email{xiao.li.xl2@roche.com}
#'
#' @export





regionalization_crosscorrelation_geostat <- function(spp_df, from_type, to_type, scale, myplot){

  L = max(max(spp_df$x_coord), max(spp_df$y_coord))
  radius = scale/L


  spp_tumor = spp_df %>% filter(cell_class == to_type) %>% select(x, y)
  spp_immune = spp_df %>% filter(cell_class == from_type) %>% select(x, y)

  # r <- raster(ncol=50, nrow=50, xmn = 0, xmx = 1, ymn = 0, ymx = 1)
  xg <- seq(0, max(spp_df$x), length.out = ceiling(max(spp_df$x)/(2*radius))+1)
  yg <- seq(0, max(spp_df$y), length.out = ceiling(max(spp_df$y)/(2*radius))+1)

  grid <- expand.grid(xg, yg)


  binxy_TC <- data.frame(x=findInterval(spp_tumor$x, xg, rightmost.closed = TRUE),
                         y=findInterval(spp_tumor$y, yg, rightmost.closed = TRUE))

  binxy_IC <- data.frame(x=findInterval(spp_immune$x, xg, rightmost.closed = TRUE),
                         y=findInterval(spp_immune$y, yg, rightmost.closed = TRUE))

  results_TC <- table(binxy_TC)
  results_IC <- table(binxy_IC)


  d_TC <- as.data.frame.table(results_TC)
  d_IC <- as.data.frame.table(results_IC)

  d_TC$x = as.numeric(as.character(d_TC$x))
  d_IC$x = as.numeric(as.character(d_IC$x))
  d_TC$y = as.numeric(as.character(d_TC$y))
  d_IC$y = as.numeric(as.character(d_IC$y))



  xx <- xg[-length(xg)] + 0.5*diff(xg)
  d_TC$xg <- xx[d_TC$x]
  d_IC$xg <- xx[d_IC$x]
  yy <- yg[-length(yg)] + 0.5*diff(yg)
  d_TC$yg <- yy[d_TC$y]
  d_IC$yg <- yy[d_IC$y]


  ### make a background grid space

  bk_grid = cbind(expand.grid(1:length(xx), 1:length(yy)), expand.grid(xx, yy))
  names(bk_grid) = c('x', 'y', 'xg', 'yg')


  d_TC = suppressMessages(left_join(bk_grid, d_TC))
  d_IC = suppressMessages(left_join(bk_grid, d_IC))

  
  ROI = !(is.na(d_TC$Freq) | is.na(d_IC$Freq) | (d_TC$Freq == 0 & d_IC$Freq == 0))
  
  d_TC$p = d_TC$Freq/nrow(spp_tumor)
  d_IC$p = d_IC$Freq/nrow(spp_immune)
  
  
  d_TC_poly = d_TC %>% filter(ROI)
  d_IC_poly = d_IC %>% filter(ROI)
  
  
  d_comb_poly = suppressMessages(left_join(d_TC_poly %>% rename(p_TC = p), 
                                           d_IC_poly %>% rename(p_IC = p),
                                           by = c("xg", "yg")) )
  
  sp::coordinates(d_comb_poly) = ~ xg + yg
  
  
  ########################################################################################################
  ############################        from from_type to to_type           ################################
  ########################################################################################################
  
  v_IC_TC = gstat::variogram(p_IC ~ p_TC, d_comb_poly)
  if (is.null(v_IC_TC)){
    v_IC_TC = gstat::variogram(p_IC ~ p_TC, d_comb_poly, cutoff = max(dist(c(d_comb_poly$xg, d_comb_poly$yg))))
    varg_IC_TC = gstat::fit.variogram(v_IC_TC, gstat::vgm("Mat"), fit.kappa = TRUE)
  } else {
    varg_IC_TC = gstat::fit.variogram(v_IC_TC, gstat::vgm("Mat"), fit.kappa = TRUE)
  }
  
  if (myplot) {
    plot_variogram(v_IC_TC, varg_IC_TC , paste0("Variogram fitting for gridded ", from_type, " and ", to_type))
  }  
  
  
  sill_IC_TC = varg_IC_TC$psill[2] + varg_IC_TC$psill[1]
  range_IC_TC = varg_IC_TC$range[2]
  kappa_IC_TC = varg_IC_TC$kappa[2]
  

  
  ########################################################################################################
  ############################        from to_type to from_type           ################################
  ########################################################################################################
  
  v_TC_IC = gstat::variogram(p_TC ~ p_IC, d_comb_poly)
  if (is.null(v_TC_IC)){
    v_TC_IC = gstat::variogram(p_TC ~ p_IC, d_comb_poly, cutoff = max(dist(c(d_comb_poly$xg, d_comb_poly$yg))))
    varg_TC_IC = gstat::fit.variogram(v_TC_IC, gstat::vgm("Mat"), fit.kappa = TRUE)
  } else {
    varg_TC_IC = gstat::fit.variogram(v_TC_IC, gstat::vgm("Mat"), fit.kappa = TRUE)
  }
  
  if (myplot) {
    plot_variogram(v_TC_IC, varg_TC_IC , paste0("Variogram fitting for gridded ", to_type, " and ", from_type))
  }  
  
  
  sill_TC_IC = varg_TC_IC$psill[2] + varg_TC_IC$psill[1]
  range_TC_IC = varg_TC_IC$range[2]
  kappa_TC_IC = varg_TC_IC$kappa[2]
  
  
  

  return(list("sill_IC_TC" = sill_IC_TC, "range_IC_TC" = range_IC_TC, "kappa_IC_TC" = kappa_IC_TC,
              "sill_TC_IC" = sill_TC_IC, "range_TC_IC" = range_TC_IC, "kappa_TC_IC" = kappa_TC_IC
             
  ))
}

