#' Semi-variogram for geostatistics data
#'
#' This is the core function that calculates the features for semi-variogram, using geostatistics methods
#'
#' @param spp_df The input spatial data.frame, need to have 3 columns: "x_coord", "y_coord" and "cell_class"
#' @param from_type The cell type one wants to use as the "from" cell type
#' @param to_type The cell type one wants to use as the "to" cell type
#' @param scale The spatial range that user wants to investigate
#' @param myplot Whether to plot the results, if available, by default it is set as FALSE
#'
#' @return This core function returns the features for semi-variogram, using geostatistics methods
#'
#' @author Xiao Li, \email{li.xiao@gene.com}
#'
#' @export





regionalization_geostat <- function(spp_df, from_type, to_type, scale,
                                    myplot){

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


  d_TC$ROI = (!(is.na(d_TC$Freq) | d_TC$Freq == 0)) | (!(is.na(d_IC$Freq) | d_IC$Freq == 0))
  d_IC$ROI = (!(is.na(d_IC$Freq) | d_IC$Freq == 0)) | (!(is.na(d_TC$Freq) | d_TC$Freq == 0))


  d_TC$p = d_TC$Freq/nrow(spp_tumor)
  d_IC$p = d_IC$Freq/nrow(spp_immune)


  d_TC_poly = d_TC %>% filter(ROI)
  d_IC_poly = d_IC %>% filter(ROI)


  ### convert the format into spatialdataframe
  coordinates(d_TC_poly) = ~ xg + yg
  coordinates(d_IC_poly) = ~ xg + yg

  if (any(is.na(d_TC_poly$p))) {
    d_TC_poly$p[is.na(d_TC_poly$p)] = 0
  }

  if (any(is.na(d_IC_poly$p))) {
    d_IC_poly$p[is.na(d_IC_poly$p)] = 0
  }

  
  
  ### Semi-variogram
  v_tumor = variogram(p~1, d_TC_poly)
  if (is.null(v_tumor)){
    v_tumor = variogram(p~1, d_TC_poly, cutoff = max(dist(c(d_TC_poly$xg, d_TC_poly$yg))))
    varg_tumor = fit.variogram(v_tumor, vgm("Mat"), fit.kappa = TRUE)
  } else {
    varg_tumor = fit.variogram(v_tumor, vgm("Mat"), fit.kappa = TRUE)
  }

  if (myplot) {
    plot_variogram(v_tumor, varg_tumor , paste0("Variogram fitting for gridded ", to_type, " proportion"))
  }  
  


  v_immune = variogram(p~1, d_IC_poly)
  if (is.null(v_immune)){
    v_immune = variogram(p~1, d_IC_poly, cutoff = max(dist(c(d_IC_poly$xg, d_IC_poly$yg))))
    varg_immune = fit.variogram(v_immune, vgm("Mat"), fit.kappa = TRUE)
  } else {
    varg_immune = fit.variogram(v_immune, vgm("Mat"), fit.kappa = TRUE)
  }

  if (myplot) {
    plot_variogram(v_immune, varg_immune , paste0("Variogram fitting for gridded ", from_type, " proportion"))
  }  
  
  
  

  sill_tumor = varg_tumor$psill[2] + varg_tumor$psill[1]
  sill_immune = varg_immune$psill[2] + varg_immune$psill[1]
  range_tumor = varg_tumor$range[2]
  range_immune = varg_immune$range[2]
  kappa_tumor = varg_tumor$kappa[2]
  kappa_immune = varg_immune$kappa[2]


  return(list("sill_tumor" = sill_tumor, "sill_immune" = sill_immune, "range_tumor" = range_tumor, "range_immune" = range_immune,
              "kappa_tumor" = kappa_tumor, "kappa_immune" = kappa_immune))
}
