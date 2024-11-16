#' Main function to generate spatial features using areal data
#'
#' This is a function that calculates the features for lattice data, using spatial lattice process.
#' Methods included in this function include global and local autocorrelation indices, e.g. Moran's I, Geary's C, Getis-Ord statistics.
#' Co-localization indices are also calculated, such as Morisita-Horn index, Sorensen index, and Jaccard index
#' 
#' 
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
#' @return This function returns the features for areal data (denote from_type as type 1 and to_type as type 2), using spatial lattice process methods:
#' \item{BC}{Bhattacharyya coefficient of type 1&2}
#' \item{MH_index}{Morisita-Horn index of type 1&2}
#' \item{JaccardJ}{Jaccard index of type 1&2}                                 
#' \item{SorensenL}{Sorensen index of type 1&2}
#' \item{Moran_I_tumor}{global Moran's I of type 1}
#' \item{Moran_I_immune}{global Moran's I of type 2}
#' \item{moran_I_Bivariate}{global bivariate Moran's I of type 1&2}
#' \item{geary_TC}{global Geary's C of type 1}
#' \item{geary_IC}{global Geary's C of type 2}
#' \item{moran_HL_TC}{percentage of significant high-low local Moran's I of type 1}
#' \item{moran_HH_TC}{percentage of significant high-high local Moran's I of type 1}
#' \item{moran_LH_TC}{percentage of significant low-high local Moran's I of type 1}
#' \item{moran_LL_TC}{percentage of significant low-low local Moran's I of type 1}
#' \item{geary_HH_TC}{percentage of significant high-high local Geary's C of type 1}
#' \item{geary_LL_TC}{percentage of significant low-low local Geary's C of type 1}
#' \item{moran_HL_IC}{percentage of significant high-low local Moran's I of type 2}
#' \item{moran_HH_IC}{percentage of significant high-high local Moran's I of type 2}
#' \item{moran_LH_IC}{percentage of significant low-high local Moran's I of type 2}
#' \item{moran_LL_IC}{percentage of significant low-low local Moran's I of type 2}
#' \item{geary_HH_IC}{percentage of significant high-high local Geary's C of type 2}
#' \item{geary_LL_IC}{percentage of significant low-low local Geary's C of type 2}
#' \item{GetisOrd_HS}{percentage of significant region being hotspot for both type 1 and type 2}
#' \item{GetisOrd_CS}{percentage of significant region being coldspot for both type 1 and type 2}
#' \item{GetisOrd_CS_IC_HS_TC}{percentage of significant region being hotspot for type 1 and coldspot for type 2}
#' \item{GetisOrd_CS_TC_HS_IC}{percentage of significant region being coldspot for type 1 and hotspot for type 2}
#' \item{GetisOrd_HS_IC}{percentage of significant region being hotspot for type 2}
#' \item{GetisOrd_CS_IC}{percentage of significant region being coldspot for type 2}
#' \item{GetisOrd_HS_TC}{percentage of significant region being hotspot for type 1}
#' \item{GetisOrd_CS_TC}{percentage of significant region being coldspot for type 1}
#' \item{GetisOrd_S_intra_cancer}{GetisOrd_HS/(GetisOrd_HS + GetisOrd_CS_IC_HS_TC + 0.000001), 0.000001 added to avoid denominator being 0}
#' \item{GetisOrd_S_intra_immune}{GetisOrd_HS/(GetisOrd_HS + GetisOrd_CS_TC_HS_IC + 0.000001), 0.000001 added to avoid denominator being 0}
#' \item{Lee_L}{global Lee's L statistics of type 1 and type 2}
#' \item{Lee_HL_TC_IC}{percentage of significant high-low local Lee's L from type 1 to type 2}
#' \item{Lee_HH_TC_IC}{percentage of significant high-high local Lee's L from type 1 to type 2}
#' \item{Lee_LH_TC_IC}{percentage of significant low-high local Lee's L from type 1 to type 2}
#' \item{Lee_LL_TC_IC}{percentage of significant low-low local Lee's L from type 1 to type 2}
#' \item{Lee_HL_IC_TC}{percentage of significant high-low local Lee's L from type 2 to type 1}
#' \item{Lee_HH_IC_TC}{percentage of significant high-high local Lee's L from type 2 to type 1}
#' \item{Lee_LH_IC_TC}{percentage of significant low-high local Lee's L from type 2 to type 1}
#' \item{Lee_LL_IC_TC}{percentage of significant low-low local Lee's L from type 2 to type 1}
#' 
#' @import magrittr
#' @import dplyr
#' @import ggplot2
#' @import spdep
#' @import sp 
#'
#' @author Xiao Li, \email{xiao.li.xl2@roche.com}
#'
#' @export





Areal_data <- function(path = "/Users/lix233/Haystack/5862_cell_centers/",
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

  
  reg_res = regionalization_lattice(spp_df = spp_df, from_type = from_type, to_type = to_type, scale = scale, myplot = myplot)
  
  reg_res_corss = regionalization_crosscorrelation_lattice(spp_df = spp_df, from_type = from_type, to_type = to_type, scale = scale, myplot = myplot)
  
  cat("Areal features are calculated. \n")

  BC = reg_res$BC
  MH_index = reg_res$MH_index
  JaccardJ = reg_res$JaccardJ
  SorensenL = reg_res$SorensenL
  moran_TC = reg_res$Moran_I_tumor
  moran_IC = reg_res$Moran_I_immune
  moran_BV = reg_res$moran_I_Bivariate
  geary_TC = reg_res$geary_TC
  geary_IC = reg_res$geary_IC
  moran_HL_TC = reg_res$moran_HL_TC
  moran_HH_TC = reg_res$moran_HH_TC
  moran_LH_TC = reg_res$moran_LH_TC
  moran_LL_TC = reg_res$moran_LL_TC
  geary_HH_TC = reg_res$geary_HH_TC
  geary_LL_TC = reg_res$geary_LL_TC
  moran_HL_IC = reg_res$moran_HL_IC
  moran_HH_IC = reg_res$moran_HH_IC
  moran_LH_IC = reg_res$moran_LH_IC
  moran_LL_IC = reg_res$moran_LL_IC
  geary_HH_IC = reg_res$geary_HH_IC
  geary_LL_IC = reg_res$geary_LL_IC
  GetisOrd_HS = reg_res$GetisOrd_HS
  GetisOrd_CS = reg_res$GetisOrd_CS
  GetisOrd_CS_IC_HS_TC = reg_res$GetisOrd_CS_IC_HS_TC
  GetisOrd_CS_TC_HS_IC = reg_res$GetisOrd_CS_TC_HS_IC
  GetisOrd_HS_IC = reg_res$GetisOrd_HS_IC
  GetisOrd_CS_IC = reg_res$GetisOrd_CS_IC
  GetisOrd_HS_TC = reg_res$GetisOrd_HS_TC
  GetisOrd_CS_TC = reg_res$GetisOrd_CS_TC
  GetisOrd_S_intra_cancer = reg_res$GetisOrd_S_intra_cancer
  GetisOrd_S_intra_immune = reg_res$GetisOrd_S_intra_immune
  
  Lee_L = reg_res_corss$Lee_L
  Lee_HL_TC_IC = reg_res_corss$Lee_HL_TC_IC
  Lee_HH_TC_IC = reg_res_corss$Lee_HH_TC_IC
  Lee_LH_TC_IC = reg_res_corss$Lee_LH_TC_IC
  Lee_LL_TC_IC = reg_res_corss$Lee_LL_TC_IC
  Lee_HL_IC_TC = reg_res_corss$Lee_HL_IC_TC
  Lee_HH_IC_TC = reg_res_corss$Lee_HH_IC_TC
  Lee_LH_IC_TC = reg_res_corss$Lee_LH_IC_TC
  Lee_LL_IC_TC = reg_res_corss$Lee_LL_IC_TC
    



  return(list("BC" = BC, "MH_index" = MH_index, "JaccardJ" = JaccardJ, "SorensenL" = SorensenL,
              "Moran_I_tumor" = moran_TC, "Moran_I_immune" = moran_IC, "moran_I_Bivariate" = moran_BV,
              "geary_TC" = geary_TC, "geary_IC" = geary_IC,
              "moran_HL_TC" = moran_HL_TC, "moran_HH_TC" = moran_HH_TC, "moran_LH_TC" = moran_LH_TC, "moran_LL_TC" = moran_LL_TC,
              "geary_HH_TC" = geary_HH_TC, "geary_LL_TC" = geary_LL_TC,
              "moran_HL_IC" = moran_HL_IC, "moran_HH_IC" = moran_HH_IC, "moran_LH_IC" = moran_LH_IC, "moran_LL_IC" = moran_LL_IC,
              "geary_HH_IC" = geary_HH_IC, "geary_LL_IC" = geary_LL_IC,
              "GetisOrd_HS" = GetisOrd_HS, "GetisOrd_CS" = GetisOrd_CS,
              "GetisOrd_CS_IC_HS_TC" = GetisOrd_CS_IC_HS_TC, "GetisOrd_CS_TC_HS_IC" = GetisOrd_CS_TC_HS_IC,
              "GetisOrd_HS_IC" = GetisOrd_HS_IC, "GetisOrd_CS_IC" = GetisOrd_CS_IC,
              "GetisOrd_HS_TC" = GetisOrd_HS_TC, "GetisOrd_CS_TC" = GetisOrd_CS_TC,
              "GetisOrd_S_intra_cancer" = GetisOrd_S_intra_cancer, "GetisOrd_S_intra_immune" = GetisOrd_S_intra_immune,
              "Lee_L" = Lee_L,
              "Lee_HL_TC_IC" = Lee_HL_TC_IC, "Lee_HH_TC_IC" = Lee_HH_TC_IC, "Lee_LH_TC_IC" = Lee_LH_TC_IC, "Lee_LL_TC_IC" = Lee_LL_TC_IC,
              "Lee_HL_IC_TC" = Lee_HL_IC_TC, "Lee_HH_IC_TC" = Lee_HH_IC_TC, "Lee_LH_IC_TC" = Lee_LH_IC_TC, "Lee_LL_IC_TC" = Lee_LL_IC_TC
              ))
}
