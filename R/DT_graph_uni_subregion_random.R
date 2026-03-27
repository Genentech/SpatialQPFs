#' The main function that calculate the Delaunay_triangulation graph features for a single cell type
#'
#' This is a function that calculate the Delaunay_triangulation graph features for a single cell type, based on random FOVs. The FOV size and number of FOVs are decided by decided by the user.
#'
#' @param path The path for the directory that contains the data
#' @param file The file name in the directory
#' @param cell_type The cell type that graph features are to be calculated
#' @param cell_class_var The column name in each file that indicates the cell class variable
#' @param x_var The column name in each file that indicates the cell location x variable
#' @param y_var The column name in each file that indicates the cell location y variable
#' @param scale The threshold to trim the graph edges. This is especially useful when the edge connects cells that are from different tissue region, e.g. different tumor nests. If user wants to retain all the edges, set it to a very large number.
#' @param side_length The side length of squared FOVs
#' @param num_FOV Number of FOVs to choose
#' @param set_seed set seed for reproducibility
#' @param myplot Whether to plot the results, if available, by default it is set as FALSE
#' 
#'
#' @return This function returns the median value of each each feature across the FOVs:
#' \item{mean_x_weight}{average edge length of the graph}
#' \item{median_x_weight}{median edge length of the graph}
#' \item{sd_x_weight}{standard deviation of edge length of the graph}
#' \item{iqr_x_weight}{IQR of edge length of the graph}
#' \item{skewness_x_weight}{skewness of edge length of the graph}
#' \item{kurtosis_x_weight}{kurtosis of edge length of the graph}
#' \item{min_x_weight}{minimum edge length of the graph}
#' \item{max_x_weight}{maximum edge length of the graph}
#' \item{range_x_weight}{range of edge length of the graph}
#' \item{Q10_weight}{10\% quantile edge length of the graph}
#' \item{Q20_weight}{20\% quantile edge length of the graph}
#' \item{Q25_weight}{25\% quantile edge length of the graph}
#' \item{Q30_weight}{30\% quantile edge length of the graph}
#' \item{Q40_weight}{40\% quantile edge length of the graph}
#' \item{Q60_weight}{60\% quantile edge length of the graph}
#' \item{Q70_weight}{70\% quantile edge length of the graph}
#' \item{Q75_weight}{75\% quantile edge length of the graph}
#' \item{Q80_weight}{80\% quantile edge length of the graph}
#' \item{Q90_weight}{90\% quantile edge length of the graph}
#' \item{theil_index_weight}{theil index edge length of the graph}
#' \item{gini_coeff_weight}{gini index edge length of the graph}
#' \item{mean_x_closeness}{average node closeness centrality of the graph}
#' \item{median_x_closeness}{median node closeness centrality of the graph}
#' \item{sd_x_closeness}{standard deviation of node closeness centrality of the graph}
#' \item{iqr_x_closeness}{IQR of node closeness centrality of the graph}
#' \item{skewness_x_closeness}{skewness of node closeness centrality of the graph}
#' \item{kurtosis_x_closeness}{kurtosis of node closeness centrality of the graph}
#' \item{min_x_closeness}{minimum node closeness centrality of the graph}
#' \item{max_x_closeness}{maximum node closeness centrality of the graph}
#' \item{range_x_closeness}{range of node closeness centrality of the graph}
#' \item{Q10_closeness}{10\% quantile node closeness centrality of the graph}
#' \item{Q20_closeness}{20\% quantile node closeness centrality of the graph}
#' \item{Q25_closeness}{25\% quantile node closeness centrality of the graph}
#' \item{Q30_closeness}{30\% quantile node closeness centrality of the graph}
#' \item{Q40_closeness}{40\% quantile node closeness centrality of the graph}
#' \item{Q60_closeness}{60\% quantile node closeness centrality of the graph}
#' \item{Q70_closeness}{70\% quantile node closeness centrality of the graph}
#' \item{Q75_closeness}{75\% quantile node closeness centrality of the graph}
#' \item{Q80_closeness}{80\% quantile node closeness centrality of the graph}
#' \item{Q90_closeness}{90\% quantile node closeness centrality of the graph}
#' \item{theil_index_closeness}{theil index node closeness centrality of the graph}
#' \item{gini_coeff_closeness}{gini index node closeness centrality of the graph}
#' \item{mean_x_betweenness}{average node betweenness centrality of the graph}
#' \item{median_x_betweenness}{median node betweenness centrality of the graph}
#' \item{sd_x_betweenness}{standard deviation of node betweenness centrality of the graph}
#' \item{iqr_x_betweenness}{IQR of node betweenness centrality of the graph}
#' \item{skewness_x_betweenness}{skewness of node betweenness centrality of the graph}
#' \item{kurtosis_x_betweenness}{kurtosis of node betweenness centrality of the graph}
#' \item{min_x_betweenness}{minimum node betweenness centrality of the graph}
#' \item{max_x_betweenness}{maximum node betweenness centrality of the graph}
#' \item{range_x_betweenness}{range of node betweenness centrality of the graph}
#' \item{Q10_betweenness}{10\% quantile node betweenness centrality of the graph}
#' \item{Q20_betweenness}{20\% quantile node betweenness centrality of the graph}
#' \item{Q25_betweenness}{25\% quantile node betweenness centrality of the graph}
#' \item{Q30_betweenness}{30\% quantile node betweenness centrality of the graph}
#' \item{Q40_betweenness}{40\% quantile node betweenness centrality of the graph}
#' \item{Q60_betweenness}{60\% quantile node betweenness centrality of the graph}
#' \item{Q70_betweenness}{70\% quantile node betweenness centrality of the graph}
#' \item{Q75_betweenness}{75\% quantile node betweenness centrality of the graph}
#' \item{Q80_betweenness}{80\% quantile node betweenness centrality of the graph}
#' \item{Q90_betweenness}{90\% quantile node betweenness centrality of the graph}
#' \item{theil_index_betweenness}{theil index node betweenness centrality of the graph}
#' \item{gini_coeff_betweenness}{gini index node betweenness centrality of the graph}
#' \item{mean_x_degree}{average node degree centrality of the graph}
#' \item{median_x_degree}{median node degree centrality of the graph}
#' \item{sd_x_degree}{standard deviation of node degree centrality of the graph}
#' \item{iqr_x_degree}{IQR of node degree centrality of the graph}
#' \item{skewness_x_degree}{skewness of node degree centrality of the graph}
#' \item{kurtosis_x_degree}{kurtosis of node degree centrality of the graph}
#' \item{min_x_degree}{minimum node degree centrality of the graph}
#' \item{max_x_degree}{maximum node degree centrality of the graph}
#' \item{range_x_degree}{range of node degree centrality of the graph}
#' \item{Q10_degree}{10\% quantile node degree centrality of the graph}
#' \item{Q20_degree}{20\% quantile node degree centrality of the graph}
#' \item{Q25_degree}{25\% quantile node degree centrality of the graph}
#' \item{Q30_degree}{30\% quantile node degree centrality of the graph}
#' \item{Q40_degree}{40\% quantile node degree centrality of the graph}
#' \item{Q60_degree}{60\% quantile node degree centrality of the graph}
#' \item{Q70_degree}{70\% quantile node degree centrality of the graph}
#' \item{Q75_degree}{75\% quantile node degree centrality of the graph}
#' \item{Q80_degree}{80\% quantile node degree centrality of the graph}
#' \item{Q90_degree}{90\% quantile node degree centrality of the graph}
#' \item{theil_index_degree}{theil index node degree centrality of the graph}
#' \item{gini_coeff_degree}{gini index node degree centrality of the graph}
#' \item{mean_x_triangle_area}{average Delaunay triangle area of the graph}
#' \item{median_x_triangle_area}{median Delaunay triangle area of the graph}
#' \item{sd_x_triangle_area}{standard deviation of Delaunay triangle area of the graph}
#' \item{iqr_x_triangle_area}{IQR of Delaunay triangle area of the graph}
#' \item{skewness_x_triangle_area}{skewness of Delaunay triangle area of the graph}
#' \item{kurtosis_x_triangle_area}{kurtosis of Delaunay triangle area of the graph}
#' \item{min_x_triangle_area}{minimum Delaunay triangle area of the graph}
#' \item{max_x_triangle_area}{maximum Delaunay triangle area of the graph}
#' \item{range_x_triangle_area}{range of Delaunay triangle area of the graph}
#' \item{Q10_triangle_area}{10\% quantile Delaunay triangle area of the graph}
#' \item{Q20_triangle_area}{20\% quantile Delaunay triangle area of the graph}
#' \item{Q25_triangle_area}{25\% quantile Delaunay triangle area of the graph}
#' \item{Q30_triangle_area}{30\% quantile Delaunay triangle area of the graph}
#' \item{Q40_triangle_area}{40\% quantile Delaunay triangle area of the graph}
#' \item{Q60_triangle_area}{60\% quantile Delaunay triangle area of the graph}
#' \item{Q70_triangle_area}{70\% quantile Delaunay triangle area of the graph}
#' \item{Q75_triangle_area}{75\% quantile Delaunay triangle area of the graph}
#' \item{Q80_triangle_area}{80\% quantile Delaunay triangle area of the graph}
#' \item{Q90_triangle_area}{90\% quantile Delaunay triangle area of the graph}
#' \item{theil_index_triangle_area}{theil index Delaunay triangle area of the graph}
#' \item{gini_coeff_triangle_area}{gini index Delaunay triangle area of the graph}
#' \item{mean_x_triangle_perimeter}{average Delaunay triangle perimeter of the graph}
#' \item{median_x_triangle_perimeter}{median Delaunay triangle perimeter of the graph}
#' \item{sd_x_triangle_perimeter}{standard deviation of Delaunay triangle perimeter of the graph}
#' \item{iqr_x_triangle_perimeter}{IQR of Delaunay triangle perimeter of the graph}
#' \item{skewness_x_triangle_perimeter}{skewness of Delaunay triangle perimeter of the graph}
#' \item{kurtosis_x_triangle_perimeter}{kurtosis of Delaunay triangle perimeter of the graph}
#' \item{min_x_triangle_perimeter}{minimum Delaunay triangle perimeter of the graph}
#' \item{max_x_triangle_perimeter}{maximum Delaunay triangle perimeter of the graph}
#' \item{range_x_triangle_perimeter}{range of Delaunay triangle perimeter of the graph}
#' \item{Q10_triangle_perimeter}{10\% quantile Delaunay triangle perimeter of the graph}
#' \item{Q20_triangle_perimeter}{20\% quantile Delaunay triangle perimeter of the graph}
#' \item{Q25_triangle_perimeter}{25\% quantile Delaunay triangle perimeter of the graph}
#' \item{Q30_triangle_perimeter}{30\% quantile Delaunay triangle perimeter of the graph}
#' \item{Q40_triangle_perimeter}{40\% quantile Delaunay triangle perimeter of the graph}
#' \item{Q60_triangle_perimeter}{60\% quantile Delaunay triangle perimeter of the graph}
#' \item{Q70_triangle_perimeter}{70\% quantile Delaunay triangle perimeter of the graph}
#' \item{Q75_triangle_perimeter}{75\% quantile Delaunay triangle perimeter of the graph}
#' \item{Q80_triangle_perimeter}{80\% quantile Delaunay triangle perimeter of the graph}
#' \item{Q90_triangle_perimeter}{90\% quantile Delaunay triangle perimeter of the graph}
#' \item{theil_index_triangle_perimeter}{theil index Delaunay triangle perimeter of the graph}
#' \item{gini_coeff_triangle_perimeter}{gini index Delaunay triangle perimeter of the graph}
#'
#'
#'
#' @import magrittr
#' @import dplyr
#' @importFrom igraph graph.data.frame set_edge_attr edges closeness betweenness degree
#' @importFrom ggplot2 ggplot aes geom_segment geom_point coord_fixed theme_void theme element_blank element_text ggtitle
#' @importFrom moments skewness kurtosis
#' @importFrom ineq Theil Gini
#'
#' @author Xiao Li, \email{xiao.li.xl2@roche.com}
#'
#' @export


DT_graph_uni_subregion_random <- function(path, file, cell_type = "lymphocyte",
                                          cell_class_var = "cell_class",
                                          x_var = "X0",
                                          y_var = "X1",
                                          scale = 1000, side_length = 2000, num_FOV = 5, set_seed = 42, myplot = FALSE){
  
  set.seed(set_seed)
  spp_df = read.csv(paste0(path, file)) %>% filter(get(cell_class_var) == cell_type)
  

  spp_df = spp_df %>% select("cell_class" = all_of(cell_class_var),
                             "x_coord" = all_of(x_var),
                             "y_coord" = all_of(y_var))

  spp_df = spp_df %>% select(cell_class, x = x_coord, y = y_coord)

  spp_df$cell_class = as.factor(spp_df$cell_class)

  L = max(max(spp_df$x), max(spp_df$y))
  # radius = scale/L

  if (nrow(spp_df %>% filter(cell_class == cell_type)) >= 3){
    if (((max(spp_df$x) - min(spp_df$x)) <= side_length) | ((max(spp_df$y) - min(spp_df$y)) <= side_length)) {

      D_res = Delaunay_triangulation(spp_df, cell_type, scale, myplot = myplot)

      res_list = list(
        'mean_x_weight' = D_res$mean_x_weight,
        'median_x_weight' = D_res$median_x_weight,
        'sd_x_weight' = D_res$sd_x_weight,
        'iqr_x_weight' = D_res$iqr_x_weight,
        'skewness_x_weight' = D_res$skewness_x_weight,
        'kurtosis_x_weight' = D_res$kurtosis_x_weight,
        'min_x_weight' = D_res$min_x_weight,
        'max_x_weight' = D_res$max_x_weight,
        'range_x_weight' = D_res$range_x_weight,
        'Q10_weight' = D_res$Q10_weight,
        'Q20_weight' = D_res$Q20_weight,
        'Q25_weight' = D_res$Q25_weight,
        'Q30_weight' = D_res$Q30_weight,
        'Q40_weight' = D_res$Q40_weight,
        'Q60_weight' = D_res$Q60_weight,
        'Q70_weight' = D_res$Q70_weight,
        'Q75_weight' = D_res$Q75_weight,
        'Q80_weight' = D_res$Q80_weight,
        'Q90_weight' = D_res$Q90_weight,
        'theil_index_weight' = D_res$theil_index_weight,
        'gini_coeff_weight' = D_res$gini_coeff_weight,
        'mean_x_closeness' = D_res$mean_x_closeness,
        'median_x_closeness' = D_res$median_x_closeness,
        'sd_x_closeness' = D_res$sd_x_closeness,
        'iqr_x_closeness' = D_res$iqr_x_closeness,
        'skewness_x_closeness' = D_res$skewness_x_closeness,
        'kurtosis_x_closeness' = D_res$kurtosis_x_closeness,
        'min_x_closeness' = D_res$min_x_closeness,
        'max_x_closeness' = D_res$max_x_closeness,
        'range_x_closeness' = D_res$range_x_closeness,
        'Q10_closeness' = D_res$Q10_closeness,
        'Q20_closeness' = D_res$Q20_closeness,
        'Q25_closeness' = D_res$Q25_closeness,
        'Q30_closeness' = D_res$Q30_closeness,
        'Q40_closeness' = D_res$Q40_closeness,
        'Q60_closeness' = D_res$Q60_closeness,
        'Q70_closeness' = D_res$Q70_closeness,
        'Q75_closeness' = D_res$Q75_closeness,
        'Q80_closeness' = D_res$Q80_closeness,
        'Q90_closeness' = D_res$Q90_closeness,
        'theil_index_closeness' = D_res$theil_index_closeness,
        'gini_coeff_closeness' = D_res$gini_coeff_closeness,
        'mean_x_betweenness' = D_res$mean_x_betweenness,
        'median_x_betweenness' = D_res$median_x_betweenness,
        'sd_x_betweenness' = D_res$sd_x_betweenness,
        'iqr_x_betweenness' = D_res$iqr_x_betweenness,
        'skewness_x_betweenness' = D_res$skewness_x_betweenness,
        'kurtosis_x_betweenness' = D_res$kurtosis_x_betweenness,
        'min_x_betweenness' = D_res$min_x_betweenness,
        'max_x_betweenness' = D_res$max_x_betweenness,
        'range_x_betweenness' = D_res$range_x_betweenness,
        'Q10_betweenness' = D_res$Q10_betweenness,
        'Q20_betweenness' = D_res$Q20_betweenness,
        'Q25_betweenness' = D_res$Q25_betweenness,
        'Q30_betweenness' = D_res$Q30_betweenness,
        'Q40_betweenness' = D_res$Q40_betweenness,
        'Q60_betweenness' = D_res$Q60_betweenness,
        'Q70_betweenness' = D_res$Q70_betweenness,
        'Q75_betweenness' = D_res$Q75_betweenness,
        'Q80_betweenness' = D_res$Q80_betweenness,
        'Q90_betweenness' = D_res$Q90_betweenness,
        'theil_index_betweenness' = D_res$theil_index_betweenness,
        'gini_coeff_betweenness' = D_res$gini_coeff_betweenness,
        'mean_x_degree' = D_res$mean_x_degree,
        'median_x_degree' = D_res$median_x_degree,
        'sd_x_degree' = D_res$sd_x_degree,
        'iqr_x_degree' = D_res$iqr_x_degree,
        'skewness_x_degree' = D_res$skewness_x_degree,
        'kurtosis_x_degree' = D_res$kurtosis_x_degree,
        'min_x_degree' = D_res$min_x_degree,
        'max_x_degree' = D_res$max_x_degree,
        'range_x_degree' = D_res$range_x_degree,
        'Q10_degree' = D_res$Q10_degree,
        'Q20_degree' = D_res$Q20_degree,
        'Q25_degree' = D_res$Q25_degree,
        'Q30_degree' = D_res$Q30_degree,
        'Q40_degree' = D_res$Q40_degree,
        'Q60_degree' = D_res$Q60_degree,
        'Q70_degree' = D_res$Q70_degree,
        'Q75_degree' = D_res$Q75_degree,
        'Q80_degree' = D_res$Q80_degree,
        'Q90_degree' = D_res$Q90_degree,
        'theil_index_degree' = D_res$theil_index_degree,
        'gini_coeff_degree' = D_res$gini_coeff_degree,
        'mean_x_triangle_area' = D_res$mean_x_triangle_area,
        'median_x_triangle_area' = D_res$median_x_triangle_area,
        'sd_x_triangle_area' = D_res$sd_x_triangle_area,
        'iqr_x_triangle_area' = D_res$iqr_x_triangle_area,
        'skewness_x_triangle_area' = D_res$skewness_x_triangle_area,
        'kurtosis_x_triangle_area' = D_res$kurtosis_x_triangle_area,
        'min_x_triangle_area' = D_res$min_x_triangle_area,
        'max_x_triangle_area' = D_res$max_x_triangle_area,
        'range_x_triangle_area' = D_res$range_x_triangle_area,
        'Q10_triangle_area' = D_res$Q10_triangle_area,
        'Q20_triangle_area' = D_res$Q20_triangle_area,
        'Q25_triangle_area' = D_res$Q25_triangle_area,
        'Q30_triangle_area' = D_res$Q30_triangle_area,
        'Q40_triangle_area' = D_res$Q40_triangle_area,
        'Q60_triangle_area' = D_res$Q60_triangle_area,
        'Q70_triangle_area' = D_res$Q70_triangle_area,
        'Q75_triangle_area' = D_res$Q75_triangle_area,
        'Q80_triangle_area' = D_res$Q80_triangle_area,
        'Q90_triangle_area' = D_res$Q90_triangle_area,
        'theil_index_triangle_area' = D_res$theil_index_triangle_area,
        'gini_coeff_triangle_area' = D_res$gini_coeff_triangle_area,
        'mean_x_triangle_perimeter' = D_res$mean_x_triangle_perimeter,
        'median_x_triangle_perimeter' = D_res$median_x_triangle_perimeter,
        'sd_x_triangle_perimeter' = D_res$sd_x_triangle_perimeter,
        'iqr_x_triangle_perimeter' = D_res$iqr_x_triangle_perimeter,
        'skewness_x_triangle_perimeter' = D_res$skewness_x_triangle_perimeter,
        'kurtosis_x_triangle_perimeter' = D_res$kurtosis_x_triangle_perimeter,
        'min_x_triangle_perimeter' = D_res$min_x_triangle_perimeter,
        'max_x_triangle_perimeter' = D_res$max_x_triangle_perimeter,
        'range_x_triangle_perimeter' = D_res$range_x_triangle_perimeter,
        'Q10_triangle_perimeter' = D_res$Q10_triangle_perimeter,
        'Q20_triangle_perimeter' = D_res$Q20_triangle_perimeter,
        'Q25_triangle_perimeter' = D_res$Q25_triangle_perimeter,
        'Q30_triangle_perimeter' = D_res$Q30_triangle_perimeter,
        'Q40_triangle_perimeter' = D_res$Q40_triangle_perimeter,
        'Q60_triangle_perimeter' = D_res$Q60_triangle_perimeter,
        'Q70_triangle_perimeter' = D_res$Q70_triangle_perimeter,
        'Q75_triangle_perimeter' = D_res$Q75_triangle_perimeter,
        'Q80_triangle_perimeter' = D_res$Q80_triangle_perimeter,
        'Q90_triangle_perimeter' = D_res$Q90_triangle_perimeter,
        'theil_index_triangle_perimeter' = D_res$theil_index_triangle_perimeter,
        'gini_coeff_triangle_perimeter' = D_res$gini_coeff_triangle_perimeter
        )


    } else {
      ## random crop the matrix (num_FOV) times and median the feature vector

      D_res = data.frame()

      for (repeat_FOV in 1:num_FOV){
        print(repeat_FOV)
        x_origin = round(runif(1, min = min(spp_df$x), max = max(spp_df$x)-side_length))
        y_origin = round(runif(1, min = min(spp_df$y), max = max(spp_df$y)-side_length))

        spp_df_sub = spp_df %>% filter(cell_class == cell_type) %>%
          filter((x >= x_origin) & (x < (x_origin + side_length)) & (y >= y_origin) & (y < (y_origin + side_length)))

        c = 0
        ## need to have at least 20 cells in the region
        while(nrow(spp_df_sub) < 3 & c < 1000) {

          x_origin = round(runif(1, min = min(spp_df$x), max = max(spp_df$x)-side_length))
          y_origin = round(runif(1, min = min(spp_df$y), max = max(spp_df$y)-side_length))

          spp_df_sub = spp_df %>% filter(cell_class == cell_type) %>%
            filter((x >= x_origin) & (x < (x_origin + side_length)) & (y >= y_origin) & (y < (y_origin + side_length)))
          # count repeated times
          c = c+1
        }

        if (c == 1000){
          next
        } else {
          this_D_res = Delaunay_triangulation(spp_df_sub, cell_type, scale, myplot = myplot)
          D_res = rbind(D_res, as.data.frame(do.call(cbind, this_D_res)))
          
        }
        
      }

      if (dim(D_res)[1] > 0) {
        D_res = data.frame(t(apply(D_res, 2, function(x) median(x[which(!is.na(x))]))))
        
        
        res_list = list(
          'mean_x_weight' = D_res$mean_x_weight,
          'median_x_weight' = D_res$median_x_weight,
          'sd_x_weight' = D_res$sd_x_weight,
          'iqr_x_weight' = D_res$iqr_x_weight,
          'skewness_x_weight' = D_res$skewness_x_weight,
          'kurtosis_x_weight' = D_res$kurtosis_x_weight,
          'min_x_weight' = D_res$min_x_weight,
          'max_x_weight' = D_res$max_x_weight,
          'range_x_weight' = D_res$range_x_weight,
          'Q10_weight' = D_res$Q10_weight,
          'Q20_weight' = D_res$Q20_weight,
          'Q25_weight' = D_res$Q25_weight,
          'Q30_weight' = D_res$Q30_weight,
          'Q40_weight' = D_res$Q40_weight,
          'Q60_weight' = D_res$Q60_weight,
          'Q70_weight' = D_res$Q70_weight,
          'Q75_weight' = D_res$Q75_weight,
          'Q80_weight' = D_res$Q80_weight,
          'Q90_weight' = D_res$Q90_weight,
          'theil_index_weight' = D_res$theil_index_weight,
          'gini_coeff_weight' = D_res$gini_coeff_weight,
          'mean_x_closeness' = D_res$mean_x_closeness,
          'median_x_closeness' = D_res$median_x_closeness,
          'sd_x_closeness' = D_res$sd_x_closeness,
          'iqr_x_closeness' = D_res$iqr_x_closeness,
          'skewness_x_closeness' = D_res$skewness_x_closeness,
          'kurtosis_x_closeness' = D_res$kurtosis_x_closeness,
          'min_x_closeness' = D_res$min_x_closeness,
          'max_x_closeness' = D_res$max_x_closeness,
          'range_x_closeness' = D_res$range_x_closeness,
          'Q10_closeness' = D_res$Q10_closeness,
          'Q20_closeness' = D_res$Q20_closeness,
          'Q25_closeness' = D_res$Q25_closeness,
          'Q30_closeness' = D_res$Q30_closeness,
          'Q40_closeness' = D_res$Q40_closeness,
          'Q60_closeness' = D_res$Q60_closeness,
          'Q70_closeness' = D_res$Q70_closeness,
          'Q75_closeness' = D_res$Q75_closeness,
          'Q80_closeness' = D_res$Q80_closeness,
          'Q90_closeness' = D_res$Q90_closeness,
          'theil_index_closeness' = D_res$theil_index_closeness,
          'gini_coeff_closeness' = D_res$gini_coeff_closeness,
          'mean_x_betweenness' = D_res$mean_x_betweenness,
          'median_x_betweenness' = D_res$median_x_betweenness,
          'sd_x_betweenness' = D_res$sd_x_betweenness,
          'iqr_x_betweenness' = D_res$iqr_x_betweenness,
          'skewness_x_betweenness' = D_res$skewness_x_betweenness,
          'kurtosis_x_betweenness' = D_res$kurtosis_x_betweenness,
          'min_x_betweenness' = D_res$min_x_betweenness,
          'max_x_betweenness' = D_res$max_x_betweenness,
          'range_x_betweenness' = D_res$range_x_betweenness,
          'Q10_betweenness' = D_res$Q10_betweenness,
          'Q20_betweenness' = D_res$Q20_betweenness,
          'Q25_betweenness' = D_res$Q25_betweenness,
          'Q30_betweenness' = D_res$Q30_betweenness,
          'Q40_betweenness' = D_res$Q40_betweenness,
          'Q60_betweenness' = D_res$Q60_betweenness,
          'Q70_betweenness' = D_res$Q70_betweenness,
          'Q75_betweenness' = D_res$Q75_betweenness,
          'Q80_betweenness' = D_res$Q80_betweenness,
          'Q90_betweenness' = D_res$Q90_betweenness,
          'theil_index_betweenness' = D_res$theil_index_betweenness,
          'gini_coeff_betweenness' = D_res$gini_coeff_betweenness,
          'mean_x_degree' = D_res$mean_x_degree,
          'median_x_degree' = D_res$median_x_degree,
          'sd_x_degree' = D_res$sd_x_degree,
          'iqr_x_degree' = D_res$iqr_x_degree,
          'skewness_x_degree' = D_res$skewness_x_degree,
          'kurtosis_x_degree' = D_res$kurtosis_x_degree,
          'min_x_degree' = D_res$min_x_degree,
          'max_x_degree' = D_res$max_x_degree,
          'range_x_degree' = D_res$range_x_degree,
          'Q10_degree' = D_res$Q10_degree,
          'Q20_degree' = D_res$Q20_degree,
          'Q25_degree' = D_res$Q25_degree,
          'Q30_degree' = D_res$Q30_degree,
          'Q40_degree' = D_res$Q40_degree,
          'Q60_degree' = D_res$Q60_degree,
          'Q70_degree' = D_res$Q70_degree,
          'Q75_degree' = D_res$Q75_degree,
          'Q80_degree' = D_res$Q80_degree,
          'Q90_degree' = D_res$Q90_degree,
          'theil_index_degree' = D_res$theil_index_degree,
          'gini_coeff_degree' = D_res$gini_coeff_degree,
          'mean_x_triangle_area' = D_res$mean_x_triangle_area,
          'median_x_triangle_area' = D_res$median_x_triangle_area,
          'sd_x_triangle_area' = D_res$sd_x_triangle_area,
          'iqr_x_triangle_area' = D_res$iqr_x_triangle_area,
          'skewness_x_triangle_area' = D_res$skewness_x_triangle_area,
          'kurtosis_x_triangle_area' = D_res$kurtosis_x_triangle_area,
          'min_x_triangle_area' = D_res$min_x_triangle_area,
          'max_x_triangle_area' = D_res$max_x_triangle_area,
          'range_x_triangle_area' = D_res$range_x_triangle_area,
          'Q10_triangle_area' = D_res$Q10_triangle_area,
          'Q20_triangle_area' = D_res$Q20_triangle_area,
          'Q25_triangle_area' = D_res$Q25_triangle_area,
          'Q30_triangle_area' = D_res$Q30_triangle_area,
          'Q40_triangle_area' = D_res$Q40_triangle_area,
          'Q60_triangle_area' = D_res$Q60_triangle_area,
          'Q70_triangle_area' = D_res$Q70_triangle_area,
          'Q75_triangle_area' = D_res$Q75_triangle_area,
          'Q80_triangle_area' = D_res$Q80_triangle_area,
          'Q90_triangle_area' = D_res$Q90_triangle_area,
          'theil_index_triangle_area' = D_res$theil_index_triangle_area,
          'gini_coeff_triangle_area' = D_res$gini_coeff_triangle_area,
          'mean_x_triangle_perimeter' = D_res$mean_x_triangle_perimeter,
          'median_x_triangle_perimeter' = D_res$median_x_triangle_perimeter,
          'sd_x_triangle_perimeter' = D_res$sd_x_triangle_perimeter,
          'iqr_x_triangle_perimeter' = D_res$iqr_x_triangle_perimeter,
          'skewness_x_triangle_perimeter' = D_res$skewness_x_triangle_perimeter,
          'kurtosis_x_triangle_perimeter' = D_res$kurtosis_x_triangle_perimeter,
          'min_x_triangle_perimeter' = D_res$min_x_triangle_perimeter,
          'max_x_triangle_perimeter' = D_res$max_x_triangle_perimeter,
          'range_x_triangle_perimeter' = D_res$range_x_triangle_perimeter,
          'Q10_triangle_perimeter' = D_res$Q10_triangle_perimeter,
          'Q20_triangle_perimeter' = D_res$Q20_triangle_perimeter,
          'Q25_triangle_perimeter' = D_res$Q25_triangle_perimeter,
          'Q30_triangle_perimeter' = D_res$Q30_triangle_perimeter,
          'Q40_triangle_perimeter' = D_res$Q40_triangle_perimeter,
          'Q60_triangle_perimeter' = D_res$Q60_triangle_perimeter,
          'Q70_triangle_perimeter' = D_res$Q70_triangle_perimeter,
          'Q75_triangle_perimeter' = D_res$Q75_triangle_perimeter,
          'Q80_triangle_perimeter' = D_res$Q80_triangle_perimeter,
          'Q90_triangle_perimeter' = D_res$Q90_triangle_perimeter,
          'theil_index_triangle_perimeter' = D_res$theil_index_triangle_perimeter,
          'gini_coeff_triangle_perimeter' = D_res$gini_coeff_triangle_perimeter
        )
      } else {
        res_list = list(
          'mean_x_weight' = NA,
          'median_x_weight' = NA,
          'sd_x_weight' = NA,
          'iqr_x_weight' = NA,
          'skewness_x_weight' = NA,
          'kurtosis_x_weight' = NA,
          'min_x_weight' = NA,
          'max_x_weight' = NA,
          'range_x_weight' = NA,
          'Q10_weight' = NA,
          'Q20_weight' = NA,
          'Q25_weight' = NA,
          'Q30_weight' = NA,
          'Q40_weight' = NA,
          'Q60_weight' = NA,
          'Q70_weight' = NA,
          'Q75_weight' = NA,
          'Q80_weight' = NA,
          'Q90_weight' = NA,
          'theil_index_weight' = NA,
          'gini_coeff_weight' = NA,
          'mean_x_closeness' = NA,
          'median_x_closeness' = NA,
          'sd_x_closeness' = NA,
          'iqr_x_closeness' = NA,
          'skewness_x_closeness' = NA,
          'kurtosis_x_closeness' = NA,
          'min_x_closeness' = NA,
          'max_x_closeness' = NA,
          'range_x_closeness' = NA,
          'Q10_closeness' = NA,
          'Q20_closeness' = NA,
          'Q25_closeness' = NA,
          'Q30_closeness' = NA,
          'Q40_closeness' = NA,
          'Q60_closeness' = NA,
          'Q70_closeness' = NA,
          'Q75_closeness' = NA,
          'Q80_closeness' = NA,
          'Q90_closeness' = NA,
          'theil_index_closeness' = NA,
          'gini_coeff_closeness' = NA,
          'mean_x_betweenness' = NA,
          'median_x_betweenness' = NA,
          'sd_x_betweenness' = NA,
          'iqr_x_betweenness' = NA,
          'skewness_x_betweenness' = NA,
          'kurtosis_x_betweenness' = NA,
          'min_x_betweenness' = NA,
          'max_x_betweenness' = NA,
          'range_x_betweenness' = NA,
          'Q10_betweenness' = NA,
          'Q20_betweenness' = NA,
          'Q25_betweenness' = NA,
          'Q30_betweenness' = NA,
          'Q40_betweenness' = NA,
          'Q60_betweenness' = NA,
          'Q70_betweenness' = NA,
          'Q75_betweenness' = NA,
          'Q80_betweenness' = NA,
          'Q90_betweenness' = NA,
          'theil_index_betweenness' = NA,
          'gini_coeff_betweenness' = NA,
          'mean_x_degree' = NA,
          'median_x_degree' = NA,
          'sd_x_degree' = NA,
          'iqr_x_degree' = NA,
          'skewness_x_degree' = NA,
          'kurtosis_x_degree' = NA,
          'min_x_degree' = NA,
          'max_x_degree' = NA,
          'range_x_degree' = NA,
          'Q10_degree' = NA,
          'Q20_degree' = NA,
          'Q25_degree' = NA,
          'Q30_degree' = NA,
          'Q40_degree' = NA,
          'Q60_degree' = NA,
          'Q70_degree' = NA,
          'Q75_degree' = NA,
          'Q80_degree' = NA,
          'Q90_degree' = NA,
          'theil_index_degree' = NA,
          'gini_coeff_degree' = NA,
          'mean_x_triangle_area' = NA,
          'median_x_triangle_area' = NA,
          'sd_x_triangle_area' = NA,
          'iqr_x_triangle_area' = NA,
          'skewness_x_triangle_area' = NA,
          'kurtosis_x_triangle_area' = NA,
          'min_x_triangle_area' = NA,
          'max_x_triangle_area' = NA,
          'range_x_triangle_area' = NA,
          'Q10_triangle_area' = NA,
          'Q20_triangle_area' = NA,
          'Q25_triangle_area' = NA,
          'Q30_triangle_area' = NA,
          'Q40_triangle_area' = NA,
          'Q60_triangle_area' = NA,
          'Q70_triangle_area' = NA,
          'Q75_triangle_area' = NA,
          'Q80_triangle_area' = NA,
          'Q90_triangle_area' = NA,
          'theil_index_triangle_area' = NA,
          'gini_coeff_triangle_area' = NA,
          'mean_x_triangle_perimeter' = NA,
          'median_x_triangle_perimeter' = NA,
          'sd_x_triangle_perimeter' = NA,
          'iqr_x_triangle_perimeter' = NA,
          'skewness_x_triangle_perimeter' = NA,
          'kurtosis_x_triangle_perimeter' = NA,
          'min_x_triangle_perimeter' = NA,
          'max_x_triangle_perimeter' = NA,
          'range_x_triangle_perimeter' = NA,
          'Q10_triangle_perimeter' = NA,
          'Q20_triangle_perimeter' = NA,
          'Q25_triangle_perimeter' = NA,
          'Q30_triangle_perimeter' = NA,
          'Q40_triangle_perimeter' = NA,
          'Q60_triangle_perimeter' = NA,
          'Q70_triangle_perimeter' = NA,
          'Q75_triangle_perimeter' = NA,
          'Q80_triangle_perimeter' = NA,
          'Q90_triangle_perimeter' = NA,
          'theil_index_triangle_perimeter' = NA,
          'gini_coeff_triangle_perimeter' = NA
        )
      }

    }
  } else {

    res_list = list(
      'mean_x_weight' = NA,
      'median_x_weight' = NA,
      'sd_x_weight' = NA,
      'iqr_x_weight' = NA,
      'skewness_x_weight' = NA,
      'kurtosis_x_weight' = NA,
      'min_x_weight' = NA,
      'max_x_weight' = NA,
      'range_x_weight' = NA,
      'Q10_weight' = NA,
      'Q20_weight' = NA,
      'Q25_weight' = NA,
      'Q30_weight' = NA,
      'Q40_weight' = NA,
      'Q60_weight' = NA,
      'Q70_weight' = NA,
      'Q75_weight' = NA,
      'Q80_weight' = NA,
      'Q90_weight' = NA,
      'theil_index_weight' = NA,
      'gini_coeff_weight' = NA,
      'mean_x_closeness' = NA,
      'median_x_closeness' = NA,
      'sd_x_closeness' = NA,
      'iqr_x_closeness' = NA,
      'skewness_x_closeness' = NA,
      'kurtosis_x_closeness' = NA,
      'min_x_closeness' = NA,
      'max_x_closeness' = NA,
      'range_x_closeness' = NA,
      'Q10_closeness' = NA,
      'Q20_closeness' = NA,
      'Q25_closeness' = NA,
      'Q30_closeness' = NA,
      'Q40_closeness' = NA,
      'Q60_closeness' = NA,
      'Q70_closeness' = NA,
      'Q75_closeness' = NA,
      'Q80_closeness' = NA,
      'Q90_closeness' = NA,
      'theil_index_closeness' = NA,
      'gini_coeff_closeness' = NA,
      'mean_x_betweenness' = NA,
      'median_x_betweenness' = NA,
      'sd_x_betweenness' = NA,
      'iqr_x_betweenness' = NA,
      'skewness_x_betweenness' = NA,
      'kurtosis_x_betweenness' = NA,
      'min_x_betweenness' = NA,
      'max_x_betweenness' = NA,
      'range_x_betweenness' = NA,
      'Q10_betweenness' = NA,
      'Q20_betweenness' = NA,
      'Q25_betweenness' = NA,
      'Q30_betweenness' = NA,
      'Q40_betweenness' = NA,
      'Q60_betweenness' = NA,
      'Q70_betweenness' = NA,
      'Q75_betweenness' = NA,
      'Q80_betweenness' = NA,
      'Q90_betweenness' = NA,
      'theil_index_betweenness' = NA,
      'gini_coeff_betweenness' = NA,
      'mean_x_degree' = NA,
      'median_x_degree' = NA,
      'sd_x_degree' = NA,
      'iqr_x_degree' = NA,
      'skewness_x_degree' = NA,
      'kurtosis_x_degree' = NA,
      'min_x_degree' = NA,
      'max_x_degree' = NA,
      'range_x_degree' = NA,
      'Q10_degree' = NA,
      'Q20_degree' = NA,
      'Q25_degree' = NA,
      'Q30_degree' = NA,
      'Q40_degree' = NA,
      'Q60_degree' = NA,
      'Q70_degree' = NA,
      'Q75_degree' = NA,
      'Q80_degree' = NA,
      'Q90_degree' = NA,
      'theil_index_degree' = NA,
      'gini_coeff_degree' = NA,
      'mean_x_triangle_area' = NA,
      'median_x_triangle_area' = NA,
      'sd_x_triangle_area' = NA,
      'iqr_x_triangle_area' = NA,
      'skewness_x_triangle_area' = NA,
      'kurtosis_x_triangle_area' = NA,
      'min_x_triangle_area' = NA,
      'max_x_triangle_area' = NA,
      'range_x_triangle_area' = NA,
      'Q10_triangle_area' = NA,
      'Q20_triangle_area' = NA,
      'Q25_triangle_area' = NA,
      'Q30_triangle_area' = NA,
      'Q40_triangle_area' = NA,
      'Q60_triangle_area' = NA,
      'Q70_triangle_area' = NA,
      'Q75_triangle_area' = NA,
      'Q80_triangle_area' = NA,
      'Q90_triangle_area' = NA,
      'theil_index_triangle_area' = NA,
      'gini_coeff_triangle_area' = NA,
      'mean_x_triangle_perimeter' = NA,
      'median_x_triangle_perimeter' = NA,
      'sd_x_triangle_perimeter' = NA,
      'iqr_x_triangle_perimeter' = NA,
      'skewness_x_triangle_perimeter' = NA,
      'kurtosis_x_triangle_perimeter' = NA,
      'min_x_triangle_perimeter' = NA,
      'max_x_triangle_perimeter' = NA,
      'range_x_triangle_perimeter' = NA,
      'Q10_triangle_perimeter' = NA,
      'Q20_triangle_perimeter' = NA,
      'Q25_triangle_perimeter' = NA,
      'Q30_triangle_perimeter' = NA,
      'Q40_triangle_perimeter' = NA,
      'Q60_triangle_perimeter' = NA,
      'Q70_triangle_perimeter' = NA,
      'Q75_triangle_perimeter' = NA,
      'Q80_triangle_perimeter' = NA,
      'Q90_triangle_perimeter' = NA,
      'theil_index_triangle_perimeter' = NA,
      'gini_coeff_triangle_perimeter' = NA
    )
  }

  return(res_list)
}
