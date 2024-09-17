#' Main function to generate spatial features using geostatistics data
#'
#' This is a function that calculates the features for point referenced data, using geostatistics methods. 
#' Semi-variogram and crossvariogram methods are implemented. 
#'
#' @param path The path for the directory that contains the data
#' @param file The file name in the directory
#' @param cell_class_var The column name in each file that indicates the cell class variable
#' @param x_var The column name in each file that indicates the cell location x variable
#' @param y_var The column name in each file that indicates the cell location y variable
#' @param from_type The cell type one wants to use as the "from" cell type
#' @param to_type The cell type one wants to use as the "to" cell type
#' @param scale The spatial range that user wants to investigate
#' @param myplot Whether to plot the results, if available, by default it is set as FALSE
#'
#' @return This function returns the features for point referenced data (denote from_type as type 1 and to_type as type 2), using geostatistics methods
#' \item{sill_tumor}{ sill of semi-variogram fitting of type 1, assuming Matérn variogram model}
#' \item{sill_immune}{ sill of semi-variogram fitting of type 2, assuming Matérn variogram model}
#' \item{range_tumor}{ range of semi-variogram fitting of type 1, assuming Matérn variogram model}
#' \item{range_immune}{ range of semi-variogram fitting of type 2, assuming Matérn variogram model}
#' \item{kappa_tumor}{ kappa of semi-variogram fitting of type 1, assuming Matérn variogram model}
#' \item{kappa_immune}{ kappa of semi-variogram fitting of type 2, assuming Matérn variogram model}
#' \item{sill_IC_TC}{ sill of cross-variogram fitting from type 2 to type 1, assuming Matérn variogram model}
#' \item{range_IC_TC}{ range of cross-variogram fitting from type 2 to type 1, assuming Matérn variogram model}
#' \item{kappa_IC_TC}{ kappa of cross-variogram fitting from type 2 to type 1, assuming Matérn variogram model}
#' \item{sill_TC_IC}{ sill of cross-variogram fitting from type 1 to type 2, assuming Matérn variogram model}
#' \item{range_TC_IC}{ range of cross-variogram fitting from type 1 to type 2, assuming Matérn variogram model}
#' \item{kappa_TC_IC}{ kappa of cross-variogram fitting from type 1 to type 2, assuming Matérn variogram model}
#' \item{kappa_IK}{ kappa of indicator-variogram fitting of type 1&2, assuming Matérn variogram model}
#' \item{sill_IK}{ sill of indicator-variogram fitting of type 1&2, assuming Matérn variogram model}
#' \item{range_IK}{ range of indicator-variogram fitting of type 1&2, assuming Matérn variogram model}
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
#' @author Xiao Li, \email{xiao.li.xl2@roche.com}
#'
#' @export





Geostatistics_data <- function(path = "/Users/lix233/Haystack/5862_cell_centers/",
                               file = "f8008cf4-3fe7-4200-a2af-21c27708ae28.csv",
                               cell_class_var = "cell_class",
                               x_var = "X0",
                               y_var = "X1",
                               from_type = "Lymphocyte",
                               to_type = "Tumor",
                               scale,
                               myplot = FALSE){

  spp_df = read.csv(paste0(path, file))

  spp_df = spp_df %>% select("cell_class" = all_of(cell_class_var),
                             "x_coord" = all_of(x_var),
                             "y_coord" = all_of(y_var))

  L = max(max(spp_df$x_coord), max(spp_df$y_coord))

  spp_df = spp_df %>% mutate("x" = normalized_coords(spp_df$x_coord, L),
                             "y" = normalized_coords(spp_df$y_coord, L),
                             "cell_class" = as.factor(spp_df$cell_class))


  radius = scale/L

  
  reg_res = regionalization_geostat(spp_df = spp_df, from_type = from_type, to_type = to_type, scale = scale, myplot = myplot)
  reg_res_corss = regionalization_crosscorrelation_geostat(spp_df = spp_df, from_type = from_type, to_type = to_type, scale = scale, myplot = myplot)
 
  sill_tumor = reg_res$sill_tumor
  sill_immune = reg_res$sill_immune
  range_tumor = reg_res$range_tumor
  range_immune = reg_res$range_immune
  kappa_tumor = reg_res$kappa_tumor
  kappa_immune = reg_res$kappa_immune
  
  sill_IC_TC = reg_res_corss$sill_IC_TC
  range_IC_TC = reg_res_corss$range_IC_TC
  kappa_IC_TC = reg_res_corss$kappa_IC_TC
  sill_TC_IC = reg_res_corss$sill_TC_IC
  range_TC_IC = reg_res_corss$range_TC_IC
  kappa_TC_IC = reg_res_corss$kappa_TC_IC
  
 
  ### Geostatistics part
  ### make the data.frame to point-referenced point data format
  spp_df$ik = spp_df$cell_class == to_type
  coordinates(spp_df) = ~ x + y

  v = variogram(ik ~ 1, spp_df)
  varg = fit.variogram(v, vgm("Mat"), fit.kappa = TRUE)

  if (myplot){
    plot_variogram(v, varg, title = "Variogram fitting for cell marks")
  }
  
  sill_IK = varg$psill[2] + varg$psill[1]
  range_IK = varg$range[2]
  kappa_IK = varg$kappa[2]

  cat("Geostatistics features are calculated. \n")



  return(list("sill_tumor" = sill_tumor, "sill_immune" = sill_immune, "range_tumor" = range_tumor, "range_immune" = range_immune,
              "kappa_tumor" = kappa_tumor, "kappa_immune" = kappa_immune,
              "sill_IC_TC" = sill_IC_TC, "range_IC_TC" = range_IC_TC, "kappa_IC_TC" = kappa_IC_TC,
              "sill_TC_IC" = sill_TC_IC, "range_TC_IC" = range_TC_IC, "kappa_TC_IC" = kappa_TC_IC,
              "kappa_IK" = kappa_IK, "sill_IK" = sill_IK, "range_IK" = range_IK))
}
