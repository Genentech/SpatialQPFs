#' The main function that calculates a rich set of Delaunay triangulation graph features.
#'
#' This function calculates Delaunay triangulation graph features for distinct cell types
#' based on randomly sampled Fields of View (FOVs). It returns a comprehensive set of
#' features describing the tissue architecture at global, local, and cellular levels.
#'
#' @param path The path for the directory that contains the data
#' @param file The file name in the directory
#' @param cell_type The cell types that graph features are to be calculated
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
#' @return This function returns a data frame with the median value of each feature across the sampled FOVs:
#'   \item{--- Global Graph Features ---}{}
#'   \item{graph_assortativity}{A score from -1 to 1. Positive values indicate spatial segregation (cells prefer same-type neighbors). Negative values indicate cell mixing.}
#'   \item{num_homotypic_edges}{The total count of edges connecting cells of the same type.}
#'   \item{num_heterotypic_edges}{The total count of edges connecting cells of different types.}
#'   \item{--- Node-Level Features (Aggregated per cell type) ---}{}
#'   \item{mean_degree_[cell_type]}{Average number of direct neighbors for a given cell type.}
#'   \item{mean_lcc_[cell_type]}{Average local clustering coefficient for a given cell type, measuring neighborhood "cliquishness".}
#'   \item{--- Neighborhood Composition Features (Interaction Quantification) ---}{}
#'   \item{avg_prop_[type_A]_neighbor_for_[type_B]}{For an average cell of type B, this is the proportion of its neighbors that are of type A.}
#'   \item{--- Original Edge Weight & Infiltration Features ---}{}
#'   \item{infiltration_score_s_k}{Ratio of heterotypic (s-k) to homotypic (s-s) edge counts.}
#'   \item{mean_x_cross_weight}{Average length of edges connecting different cell types.}
#'   \item{... (and all other edge weight statistics from your original script) ...}{}
#'
#'#' @import magrittr
#' @import dplyr
#' @importFrom igraph graph_from_data_frame E V delete_edges as_data_frame degree transitivity assortativity_nominal neighborhood
#' @importFrom RTriangle triangulate pslg
#' @importFrom ggplot2 ggplot aes geom_segment geom_point coord_fixed theme_void theme element_blank element_text ggtitle
#' @importFrom moments skewness kurtosis
#' @importFrom ineq Theil Gini
#' @importFrom stats dist
#' @importFrom data.table setnames
#'
#' @author Xiao Li, \email{xiao.li.xl2@roche.com}
#' 
#' @export



DT_graph_cross_subregion_random <- function(path, file, cell_type = c("tumor", "lymphocyte"),
                                            cell_class_var = "cell_class",
                                            x_var = "X0", y_var = "X1",
                                            scale = 1000, side_length = 2000, num_FOV = 5,
                                            set_seed = 42, myplot = FALSE){
  
  set.seed(set_seed)
  spp_df <- read.csv(paste0(path, file)) %>%
    select(
      "cell_class" = all_of(cell_class_var),
      "x" = all_of(x_var),
      "y" = all_of(y_var)
    ) %>%
    mutate(cell_class = as.factor(cell_class))
  
  if (nrow(spp_df %>% filter(cell_class %in% cell_type)) < 4) {
    warning("Not enough cells of specified types in the entire file.")
    return(NULL)
  }
  
  if (((max(spp_df$x) - min(spp_df$x)) <= side_length) || ((max(spp_df$y) - min(spp_df$y)) <= side_length)) {
    # Image is small, process the whole thing once
    res_list <- Delaunay_triangulation_cross(spp_df, cell_type, scale, myplot = myplot)
    
  } else {
    # Image is large, sample multiple FOVs
    fov_results <- list()
    
    for (repeat_FOV in 1:num_FOV) {
      cat(paste("Processing FOV:", repeat_FOV, "\n"))
      spp_df_sub <- data.frame()
      c <- 0
      while(nrow(spp_df_sub) < 4 && c < 1000) {
        x_origin <- runif(1, min = min(spp_df$x), max = max(spp_df$x) - side_length)
        y_origin <- runif(1, min = min(spp_df$y), max = max(spp_df$y) - side_length)
        spp_df_sub <- spp_df %>%
          filter(
            cell_class %in% cell_type,
            x >= x_origin & x < (x_origin + side_length) &
              y >= y_origin & y < (y_origin + side_length)
          )
        c <- c + 1
      }
      
      if (c == 1000) {
        warning(paste("Could not find a valid FOV with at least 4 cells in iteration", repeat_FOV))
        next
      }
      
      this_D_res <- Delaunay_triangulation_cross(spp_df_sub, cell_type, scale, myplot = myplot)
      if (!is.null(this_D_res)) {
        fov_results[[length(fov_results) + 1]] <- this_D_res
      }
    }
    
    if (length(fov_results) > 0) {
      # Use bind_rows for robustly combining lists into a data frame
      D_res_df <- bind_rows(fov_results)
      # Aggregate by taking the median of each feature column
      res_list <- as.list(apply(D_res_df, 2, median, na.rm = TRUE))
    } else {
      warning("No valid FOVs were successfully processed.")
      return(NULL)
    }
  }
  
  return(as.data.frame(res_list))
}