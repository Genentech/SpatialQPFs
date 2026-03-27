#' A cross type Delaunay triangulation function with comprehensive feature extraction
#'
#' This function builds a Delaunay graph and calculates a rich set of features describing
#' the spatial architecture, including global, node-level, and neighborhood metrics.
#'
#' @param spp_df The input spatial data.frame, needs columns: "x", "y", "cell_class".
#' @param cell_type The cell types that graph features are to be calculated for.
#' @param scale The threshold to trim long graph edges.
#' @param myplot Whether to plot the results (default: FALSE).
#'
#' @return A named list containing a comprehensive set of graph features.
#' (See documentation for DT_graph_cross_subregion_random for a full list).
#'
#' @import magrittr
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
#' @export
#' 



Delaunay_triangulation_cross <- function(spp_df, cell_type, scale, myplot = FALSE) {
  
  df <- spp_df %>%
    filter(cell_class %in% cell_type) %>%
    distinct(x, y, .keep_all = TRUE) # Cleaner way to remove duplicates

  
  if (nrow(df) < 4) {
    # Return a named list of NAs for consistency if graph can't be built
    # This part can be expanded to match all possible output names
    return(list(mean_x_cross_weight = NA, graph_assortativity = NA))
  }
  
  df$vertex_id <- 1:nrow(df)
  df <- df %>% select(vertex_id, everything())
  
  
  # --- Original Graph Construction using RTriangle ---
  DV_graph <- RTriangle::triangulate(RTriangle::pslg(df %>% select(x, y)))
  
  # --- PLOTTING FUNCTIONALITY RESTORED HERE ---
  if(myplot){
    
    df_nodes <- data.frame(
      x = DV_graph$P[,1],
      y = DV_graph$P[,2],
      class = df$cell_class
    )
    
    triangles <- DV_graph$T
    
    df_edges <- lapply(seq_len(nrow(triangles)), function(r) {
      tri <- triangles[r,]
      edges <- rbind(
        c(tri[1], tri[2]),
        c(tri[2], tri[3]),
        c(tri[3], tri[1])
      )
      edges
    }) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    colnames(df_edges) <- c("i1", "i2")
    
    df_edges <- df_edges %>%
      mutate(x1 = DV_graph$P[i1,1], y1 = DV_graph$P[i1,2],
             x2 = DV_graph$P[i2,1], y2 = DV_graph$P[i2,2])
    
    
    fig <- ggplot2::ggplot() +
      ggplot2::geom_segment(data = df_edges, ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
                   color = "gray50") +
      ggplot2::geom_point(data = df_nodes, ggplot2::aes(x = x, y = y, color = class), size = 3) +
      ggplot2::coord_fixed() +
      ggplot2::theme_void() +
      ggplot2::theme(legend.position = "right",
            legend.title = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")) +
      ggplot2::ggtitle("Delaunay Triangulation")
    
    print(fig)
    
  }
  # --- END OF PLOTTING CODE ---
  
  
  g_df <- data.frame(A = DV_graph$E[,1], B = DV_graph$E[,2])
  g <- igraph::graph_from_data_frame(d = g_df, directed = FALSE, vertices = df)
  
  edge_weights <- sqrt((df$x[g_df$A] - df$x[g_df$B])^2 + (df$y[g_df$A] - df$y[g_df$B])^2)
  igraph::E(g)$weight <- edge_weights
  
  # Prune the network based on the scale threshold
  g_pruned <- igraph::delete_edges(g, igraph::E(g)[igraph::E(g)$weight > scale])
  
  # Exit if no edges remain after pruning
  if (igraph::gsize(g_pruned) == 0) return(NULL)
  
  # --- Start of New Feature Extraction ---
  features <- list()
  edge_df <- igraph::as_data_frame(g_pruned, what = "edges")
  edge_df$type1 <- igraph::V(g_pruned)$cell_class[as.integer(edge_df$from)]
  edge_df$type2 <- igraph::V(g_pruned)$cell_class[as.integer(edge_df$to)]
  
  # 1. GLOBAL GRAPH FEATURES
  features$num_homotypic_edges <- sum(edge_df$type1 == edge_df$type2)
  features$num_heterotypic_edges <- sum(edge_df$type1 != edge_df$type2)
  
  if (length(unique(igraph::V(g_pruned)$cell_class)) > 1) {
    features$graph_assortativity <- igraph::assortativity_nominal(g_pruned, types = as.numeric(as.factor(igraph::V(g_pruned)$cell_class)), directed = FALSE)
  } else {
    features$graph_assortativity <- NA
  }
  
  # 2. AGGREGATED NODE-LEVEL FEATURES
  igraph::V(g_pruned)$degree <- igraph::degree(g_pruned)
  igraph::V(g_pruned)$lcc <- igraph::transitivity(g_pruned, type = "local", isolates = "zero") # Local Clustering Coeff
  
  for (ct in cell_type) {
    nodes_of_type <- igraph::V(g_pruned)[igraph::V(g_pruned)$cell_class == ct]
    if (length(nodes_of_type) > 0) {
      features[[paste0("mean_degree_", ct)]] <- mean(nodes_of_type$degree, na.rm = TRUE)
      features[[paste0("mean_lcc_", ct)]] <- mean(nodes_of_type$lcc, na.rm = TRUE)
    } else {
      features[[paste0("mean_degree_", ct)]] <- NA
      features[[paste0("mean_lcc_", ct)]] <- NA
    }
  }
  
  # 3. NEIGHBORHOOD COMPOSITION FEATURES (very powerful!)
  all_nodes <- igraph::V(g_pruned)
  neighbors_list <- igraph::neighborhood(g_pruned, order = 1, nodes = all_nodes)
  
  for (node_type_B in cell_type) {
    indices_B <- which(igraph::V(g_pruned)$cell_class == node_type_B)
    if (length(indices_B) == 0) next
    
    prop_list <- lapply(indices_B, function(i) {
      neighbor_nodes <- neighbors_list[[i]][-1]
      if (length(neighbor_nodes) > 0) {
        neighbor_classes <- igraph::V(g_pruned)$cell_class[neighbor_nodes]
        prop.table(table(factor(neighbor_classes, levels = cell_type)))
      }
    })
    prop_list <- prop_list[!sapply(prop_list, is.null)]
    
    if (length(prop_list) > 0) {
      avg_props <- colMeans(do.call(rbind, prop_list), na.rm = TRUE)
      for (node_type_A in cell_type) {
        feature_name <- paste0("avg_prop_", node_type_A, "_neighbor_for_", node_type_B)
        features[[feature_name]] <- avg_props[node_type_A]
      }
    }
  }
  
  # --- Original Infiltration Score (calculated more efficiently) ---
  for (k in cell_type) {
    for (s in setdiff(cell_type, k)) {
      num_sk <- sum((edge_df$type1 == s & edge_df$type2 == k) | (edge_df$type1 == k & edge_df$type2 == s))
      num_ss <- sum(edge_df$type1 == s & edge_df$type2 == s)
      score_name <- paste0("infiltration_score_", k, "_", s)
      features[[score_name]] <- ifelse(num_ss > 0, num_sk / num_ss, NA)
    }
  }
  
  # --- Original Edge Weight Distribution Features ---
  all_weights <- edge_df %>% filter(type1 != type2) %>% .$weight
  dist_summary <- list()
  if (length(all_weights) > 1) {
    dist_summary <- distribution_summary(all_weights)
    # The setnames function requires data.table, this is a base R equivalent
    names(dist_summary) <- paste0(names(dist_summary), "_cross_weight")
  }
  
  return(append(features, dist_summary))
}


