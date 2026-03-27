#' A function for calculating Delaunay_triangulation graph features for one cell type
#'
#' This is a function that calculate the graph features from Delaunay_triangulation graph for one cell type
#'
#' @param spp_df The input spatial data.frame, need to have 3 columns: "x", "y" and "cell_class"
#' @param cell_type The cell type that graph features are to be calculated
#' @param scale The threshold to trim the graph edges
#' @param myplot Whether to plot the results, if available, by default it is set as FALSE
#' 
#' 
#' @return This function returns the features from Delaunay_triangulation graph:
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
#' @import magrittr
#' @import dplyr
#' @importFrom igraph graph.data.frame set_edge_attr edges closeness betweenness degree 
#' @importFrom plyr join.keys
#' @importFrom ggplot2 ggplot aes geom_segment geom_point coord_fixed theme_void theme element_blank element_text ggtitle
#' @importFrom moments skewness kurtosis
#' @importFrom ineq Theil Gini
#'
#' @author Xiao Li, \email{xiao.li.xl2@roche.com}
#'
#' @export


Delaunay_triangulation <- function(spp_df, cell_type, scale, myplot) {


  df <- spp_df %>% filter(cell_class == cell_type) %>% select(x, y) %>% as.matrix()

  
  if (nrow(df) > 1){
    df <- df[!duplicated(df[,c('x', 'y')]),]
  }


  if (nrow(df) >= 3) {

    DV_graph = RTriangle::triangulate(RTriangle::pslg(df))

    if(myplot){
      
      df_nodes <- data.frame(
        x = DV_graph$P[,1],
        y = DV_graph$P[,2],
        class = rep(cell_type, nrow(df))
      )
      
      triangles <- DV_graph$T
      
      df_edges <- lapply(seq_len(nrow(triangles)), function(r) {
        tri <- triangles[r,]
        # The vertex indices for this triangle
        # We'll form edges between these points:
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
        # Draw edges (use a neutral color for edges)
        ggplot2::geom_segment(data = df_edges, ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
                     color = "gray50") +
        # Draw nodes, colored by class
        ggplot2::geom_point(data = df_nodes, ggplot2::aes(x = x, y = y, color = class), size = 3) +
        # Equal aspect ratio is often nice for geometric plots
        ggplot2::coord_fixed() +
        # Remove axes and background grid
        ggplot2::theme_void() +
        # Move legend outside the plot; for example, to the right:
        ggplot2::theme(legend.position = "right",
              legend.title = ggplot2::element_blank(),
              # Center the title
              plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")) +
        ggplot2::ggtitle("Delaunay Triangulation")
      
      print(fig)
      
    }
    
    g_df <- data.frame(
      A = DV_graph$E[,1] ,
      B = DV_graph$E[,2])

    g <- igraph::graph.data.frame(d = g_df, directed = FALSE)

    D_edges = data.frame()
    for (i in 1:nrow(DV_graph$E)) {
      D_edges = rbind(D_edges, data.frame("from" = DV_graph$E[i,1], "to" = DV_graph$E[i,2], "edge_length" = dist(df[DV_graph$E[i,], c('x', 'y')])[1]))
    }

    g <- g %>% igraph::set_edge_attr("weight", value = D_edges$edge_length)

    ## prune the network
    ### Convert the data.frame to edges
    e <- apply(D_edges[D_edges$edge_length > scale, c("from", "to")], 1, paste, collapse="|")
    e <- igraph::edges(e)

    g_pruned = g-e

    closeness_cen = igraph::closeness(g_pruned, normalized = T)
    betweenness_cen = igraph::betweenness(g_pruned, normalized = T)
    degree_cen = igraph::degree(g_pruned, normalized = T)

    tmp = D_edges[D_edges$edge_length < scale, 'edge_length']
    
    weight_df <- distribution_summary(tmp) %>% setnames(paste0(names(.),"_weight"))
    
    closeness_df <- distribution_summary(closeness_cen) %>% setnames(paste0(names(.),"_closeness"))
    
    betweenness_df <- distribution_summary(betweenness_cen) %>% setnames(paste0(names(.),"_betweenness"))

    degree_df <- distribution_summary(degree_cen) %>% setnames(paste0(names(.),"_degree"))
    
    
    mean_x_weight = weight_df$mean_x_weight 
    median_x_weight = weight_df$median_x_weight
    sd_x_weight = weight_df$sd_x_weight
    iqr_x_weight = weight_df$iqr_x_weight
    skewness_x_weight = weight_df$skewness_x_weight
    kurtosis_x_weight = weight_df$kurtosis_x_weight
    min_x_weight = weight_df$min_x_weight
    max_x_weight = weight_df$max_x_weight
    range_x_weight = weight_df$range_x_weight
    Q10_weight = weight_df$Q10_weight
    Q20_weight = weight_df$Q20_weight
    Q25_weight = weight_df$Q25_weight
    Q30_weight = weight_df$Q30_weight
    Q40_weight = weight_df$Q40_weight
    Q60_weight = weight_df$Q60_weight
    Q70_weight = weight_df$Q70_weight
    Q75_weight = weight_df$Q75_weight
    Q80_weight = weight_df$Q80_weight
    Q90_weight = weight_df$Q90_weight
    theil_index_weight = weight_df$theil_index_weight
    gini_coeff_weight = weight_df$gini_coeff_weight
    
    mean_x_closeness = closeness_df$mean_x_closeness 
    median_x_closeness = closeness_df$median_x_closeness
    sd_x_closeness = closeness_df$sd_x_closeness
    iqr_x_closeness = closeness_df$iqr_x_closeness
    skewness_x_closeness = closeness_df$skewness_x_closeness
    kurtosis_x_closeness = closeness_df$kurtosis_x_closeness
    min_x_closeness = closeness_df$min_x_closeness
    max_x_closeness = closeness_df$max_x_closeness
    range_x_closeness = closeness_df$range_x_closeness
    Q10_closeness = closeness_df$Q10_closeness
    Q20_closeness = closeness_df$Q20_closeness
    Q25_closeness = closeness_df$Q25_closeness
    Q30_closeness = closeness_df$Q30_closeness
    Q40_closeness = closeness_df$Q40_closeness
    Q60_closeness = closeness_df$Q60_closeness
    Q70_closeness = closeness_df$Q70_closeness
    Q75_closeness = closeness_df$Q75_closeness
    Q80_closeness = closeness_df$Q80_closeness
    Q90_closeness = closeness_df$Q90_closeness
    theil_index_closeness = closeness_df$theil_index_closeness
    gini_coeff_closeness = closeness_df$gini_coeff_closeness
    
    mean_x_betweenness = betweenness_df$mean_x_betweenness 
    median_x_betweenness = betweenness_df$median_x_betweenness
    sd_x_betweenness = betweenness_df$sd_x_betweenness
    iqr_x_betweenness = betweenness_df$iqr_x_betweenness
    skewness_x_betweenness = betweenness_df$skewness_x_betweenness
    kurtosis_x_betweenness = betweenness_df$kurtosis_x_betweenness
    min_x_betweenness = betweenness_df$min_x_betweenness
    max_x_betweenness = betweenness_df$max_x_betweenness
    range_x_betweenness = betweenness_df$range_x_betweenness
    Q10_betweenness = betweenness_df$Q10_betweenness
    Q20_betweenness = betweenness_df$Q20_betweenness
    Q25_betweenness = betweenness_df$Q25_betweenness
    Q30_betweenness = betweenness_df$Q30_betweenness
    Q40_betweenness = betweenness_df$Q40_betweenness
    Q60_betweenness = betweenness_df$Q60_betweenness
    Q70_betweenness = betweenness_df$Q70_betweenness
    Q75_betweenness = betweenness_df$Q75_betweenness
    Q80_betweenness = betweenness_df$Q80_betweenness
    Q90_betweenness = betweenness_df$Q90_betweenness
    theil_index_betweenness = betweenness_df$theil_index_betweenness
    gini_coeff_betweenness = betweenness_df$gini_coeff_betweenness
    
    mean_x_degree = degree_df$mean_x_degree 
    median_x_degree = degree_df$median_x_degree
    sd_x_degree = degree_df$sd_x_degree
    iqr_x_degree = degree_df$iqr_x_degree
    skewness_x_degree = degree_df$skewness_x_degree
    kurtosis_x_degree = degree_df$kurtosis_x_degree
    min_x_degree = degree_df$min_x_degree
    max_x_degree = degree_df$max_x_degree
    range_x_degree = degree_df$range_x_degree
    Q10_degree = degree_df$Q10_degree
    Q20_degree = degree_df$Q20_degree
    Q25_degree = degree_df$Q25_degree
    Q30_degree = degree_df$Q30_degree
    Q40_degree = degree_df$Q40_degree
    Q60_degree = degree_df$Q60_degree
    Q70_degree = degree_df$Q70_degree
    Q75_degree = degree_df$Q75_degree
    Q80_degree = degree_df$Q80_degree
    Q90_degree = degree_df$Q90_degree
    theil_index_degree = degree_df$theil_index_degree
    gini_coeff_degree = degree_df$gini_coeff_degree
    
    
    
    ## statistical descriptors

    rm_vertics = D_edges[D_edges$edge_length > scale, c("from", "to")]


    df1 = DV_graph$T[, 1:2]
    df2 = DV_graph$T[, 2:3]
    df3 = DV_graph$T[, c(1,3)]

    df1 = as.data.frame(df1); names(df1) = c("from", "to")
    df2 = as.data.frame(df2); names(df2) = c("from", "to")
    df3 = as.data.frame(df3); names(df3) = c("from", "to")


    re_df = rbind(as.matrix(rm_vertics), cbind(as.matrix(rm_vertics)[,2], as.matrix(rm_vertics)[,1]))

    ## we want to record the row numbers in df1, df2, df3 in order to remove them in DV_graph$T
    rm_trig = c(with(plyr::join.keys(df1, as.data.frame(re_df)), which(x %in% y)),
                with(plyr::join.keys(df2, as.data.frame(re_df)), which(x %in% y)),
                with(plyr::join.keys(df3, as.data.frame(re_df)), which(x %in% y)))


    if (length(rm_trig) > 0){
      DV_T_prune = as.data.frame(DV_graph$T)[-rm_trig, ]
    } else {
      DV_T_prune = as.data.frame(DV_graph$T)
    }


    if (nrow(DV_T_prune) >0){
      ## area of triangles
      Area_tri = c()

      for (k in 1:nrow(DV_T_prune)) {

        coords_Tri = df[as.numeric(DV_T_prune[k,]), c('x', 'y')]

        wk_mat = cbind(as.matrix(coords_Tri), rep(1, 3))

        Area_tri = c(Area_tri, 0.5*abs(base::det(wk_mat)))

      }

      
      Area_df <- distribution_summary(Area_tri) %>% setnames(paste0(names(.),"_triangle_area"))
      

      ## perimeter of triangles
      Peri_tri = c()

      for (k in 1:nrow(DV_T_prune)) {

        coords_Tri = df[as.numeric(DV_T_prune[k,]), c('x', 'y')]

        Peri_tri = c(Peri_tri, sum(dist(coords_Tri)))

      }

      Peri_df <- distribution_summary(Peri_tri) %>% setnames(paste0(names(.),"_triangle_perimeter"))
      
      
      mean_x_triangle_area = Area_df$mean_x_triangle_area 
      median_x_triangle_area = Area_df$median_x_triangle_area
      sd_x_triangle_area = Area_df$sd_x_triangle_area
      iqr_x_triangle_area = Area_df$iqr_x_triangle_area
      skewness_x_triangle_area = Area_df$skewness_x_triangle_area
      kurtosis_x_triangle_area = Area_df$kurtosis_x_triangle_area
      min_x_triangle_area = Area_df$min_x_triangle_area
      max_x_triangle_area = Area_df$max_x_triangle_area
      range_x_triangle_area = Area_df$range_x_triangle_area
      Q10_triangle_area = Area_df$Q10_triangle_area
      Q20_triangle_area = Area_df$Q20_triangle_area
      Q25_triangle_area = Area_df$Q25_triangle_area
      Q30_triangle_area = Area_df$Q30_triangle_area
      Q40_triangle_area = Area_df$Q40_triangle_area
      Q60_triangle_area = Area_df$Q60_triangle_area
      Q70_triangle_area = Area_df$Q70_triangle_area
      Q75_triangle_area = Area_df$Q75_triangle_area
      Q80_triangle_area = Area_df$Q80_triangle_area
      Q90_triangle_area = Area_df$Q90_triangle_area
      theil_index_triangle_area = Area_df$theil_index_triangle_area
      gini_coeff_triangle_area = Area_df$gini_coeff_triangle_area
      
      mean_x_triangle_perimeter = Peri_df$mean_x_triangle_perimeter 
      median_x_triangle_perimeter = Peri_df$median_x_triangle_perimeter
      sd_x_triangle_perimeter = Peri_df$sd_x_triangle_perimeter
      iqr_x_triangle_perimeter = Peri_df$iqr_x_triangle_perimeter
      skewness_x_triangle_perimeter = Peri_df$skewness_x_triangle_perimeter
      kurtosis_x_triangle_perimeter = Peri_df$kurtosis_x_triangle_perimeter
      min_x_triangle_perimeter = Peri_df$min_x_triangle_perimeter
      max_x_triangle_perimeter = Peri_df$max_x_triangle_perimeter
      range_x_triangle_perimeter = Peri_df$range_x_triangle_perimeter
      Q10_triangle_perimeter = Peri_df$Q10_triangle_perimeter
      Q20_triangle_perimeter = Peri_df$Q20_triangle_perimeter
      Q25_triangle_perimeter = Peri_df$Q25_triangle_perimeter
      Q30_triangle_perimeter = Peri_df$Q30_triangle_perimeter
      Q40_triangle_perimeter = Peri_df$Q40_triangle_perimeter
      Q60_triangle_perimeter = Peri_df$Q60_triangle_perimeter
      Q70_triangle_perimeter = Peri_df$Q70_triangle_perimeter
      Q75_triangle_perimeter = Peri_df$Q75_triangle_perimeter
      Q80_triangle_perimeter = Peri_df$Q80_triangle_perimeter
      Q90_triangle_perimeter = Peri_df$Q90_triangle_perimeter
      theil_index_triangle_perimeter = Peri_df$theil_index_triangle_perimeter
      gini_coeff_triangle_perimeter = Peri_df$gini_coeff_triangle_perimeter
      
      
      
    } else {

      mean_x_triangle_area = NA 
      median_x_triangle_area = NA
      sd_x_triangle_area = NA
      iqr_x_triangle_area = NA
      skewness_x_triangle_area = NA
      kurtosis_x_triangle_area = NA
      min_x_triangle_area = NA
      max_x_triangle_area = NA
      range_x_triangle_area = NA
      Q10_triangle_area = NA
      Q20_triangle_area = NA
      Q25_triangle_area = NA
      Q30_triangle_area = NA
      Q40_triangle_area = NA
      Q60_triangle_area = NA
      Q70_triangle_area = NA
      Q75_triangle_area = NA
      Q80_triangle_area = NA
      Q90_triangle_area = NA
      theil_index_triangle_area = NA
      gini_coeff_triangle_area = NA
      
      mean_x_triangle_perimeter = NA 
      median_x_triangle_perimeter = NA
      sd_x_triangle_perimeter = NA
      iqr_x_triangle_perimeter = NA
      skewness_x_triangle_perimeter = NA
      kurtosis_x_triangle_perimeter = NA
      min_x_triangle_perimeter = NA
      max_x_triangle_perimeter = NA
      range_x_triangle_perimeter = NA
      Q10_triangle_perimeter = NA
      Q20_triangle_perimeter = NA
      Q25_triangle_perimeter = NA
      Q30_triangle_perimeter = NA
      Q40_triangle_perimeter = NA
      Q60_triangle_perimeter = NA
      Q70_triangle_perimeter = NA
      Q75_triangle_perimeter = NA
      Q80_triangle_perimeter = NA
      Q90_triangle_perimeter = NA
      theil_index_triangle_perimeter = NA
      gini_coeff_triangle_perimeter = NA

    }

  } else {
    mean_x_weight = NA 
    median_x_weight = NA
    sd_x_weight = NA
    iqr_x_weight = NA
    skewness_x_weight = NA
    kurtosis_x_weight = NA
    min_x_weight = NA
    max_x_weight = NA
    range_x_weight = NA
    Q10_weight = NA
    Q20_weight = NA
    Q25_weight = NA
    Q30_weight = NA
    Q40_weight = NA
    Q60_weight = NA
    Q70_weight = NA
    Q75_weight = NA
    Q80_weight = NA
    Q90_weight = NA
    theil_index_weight = NA
    gini_coeff_weight = NA
    
    mean_x_closeness = NA 
    median_x_closeness = NA
    sd_x_closeness = NA
    iqr_x_closeness = NA
    skewness_x_closeness = NA
    kurtosis_x_closeness = NA
    min_x_closeness = NA
    max_x_closeness = NA
    range_x_closeness = NA
    Q10_closeness = NA
    Q20_closeness = NA
    Q25_closeness = NA
    Q30_closeness = NA
    Q40_closeness = NA
    Q60_closeness = NA
    Q70_closeness = NA
    Q75_closeness = NA
    Q80_closeness = NA
    Q90_closeness = NA
    theil_index_closeness = NA
    gini_coeff_closeness = NA
    
    mean_x_betweenness = NA 
    median_x_betweenness = NA
    sd_x_betweenness = NA
    iqr_x_betweenness = NA
    skewness_x_betweenness = NA
    kurtosis_x_betweenness = NA
    min_x_betweenness = NA
    max_x_betweenness = NA
    range_x_betweenness = NA
    Q10_betweenness = NA
    Q20_betweenness = NA
    Q25_betweenness = NA
    Q30_betweenness = NA
    Q40_betweenness = NA
    Q60_betweenness = NA
    Q70_betweenness = NA
    Q75_betweenness = NA
    Q80_betweenness = NA
    Q90_betweenness = NA
    theil_index_betweenness = NA
    gini_coeff_betweenness = NA
    
    mean_x_degree = NA 
    median_x_degree = NA
    sd_x_degree = NA
    iqr_x_degree = NA
    skewness_x_degree = NA
    kurtosis_x_degree = NA
    min_x_degree = NA
    max_x_degree = NA
    range_x_degree = NA
    Q10_degree = NA
    Q20_degree = NA
    Q25_degree = NA
    Q30_degree = NA
    Q40_degree = NA
    Q60_degree = NA
    Q70_degree = NA
    Q75_degree = NA
    Q80_degree = NA
    Q90_degree = NA
    theil_index_degree = NA
    gini_coeff_degree = NA
    
    mean_x_triangle_area = NA 
    median_x_triangle_area = NA
    sd_x_triangle_area = NA
    iqr_x_triangle_area = NA
    skewness_x_triangle_area = NA
    kurtosis_x_triangle_area = NA
    min_x_triangle_area = NA
    max_x_triangle_area = NA
    range_x_triangle_area = NA
    Q10_triangle_area = NA
    Q20_triangle_area = NA
    Q25_triangle_area = NA
    Q30_triangle_area = NA
    Q40_triangle_area = NA
    Q60_triangle_area = NA
    Q70_triangle_area = NA
    Q75_triangle_area = NA
    Q80_triangle_area = NA
    Q90_triangle_area = NA
    theil_index_triangle_area = NA
    gini_coeff_triangle_area = NA
    
    mean_x_triangle_perimeter = NA 
    median_x_triangle_perimeter = NA
    sd_x_triangle_perimeter = NA
    iqr_x_triangle_perimeter = NA
    skewness_x_triangle_perimeter = NA
    kurtosis_x_triangle_perimeter = NA
    min_x_triangle_perimeter = NA
    max_x_triangle_perimeter = NA
    range_x_triangle_perimeter = NA
    Q10_triangle_perimeter = NA
    Q20_triangle_perimeter = NA
    Q25_triangle_perimeter = NA
    Q30_triangle_perimeter = NA
    Q40_triangle_perimeter = NA
    Q60_triangle_perimeter = NA
    Q70_triangle_perimeter = NA
    Q75_triangle_perimeter = NA
    Q80_triangle_perimeter = NA
    Q90_triangle_perimeter = NA
    theil_index_triangle_perimeter = NA
    gini_coeff_triangle_perimeter = NA
  }






  return(list(
    'mean_x_weight' = mean_x_weight,
    'median_x_weight' = median_x_weight,
    'sd_x_weight' = sd_x_weight,
    'iqr_x_weight' = iqr_x_weight,
    'skewness_x_weight' = skewness_x_weight,
    'kurtosis_x_weight' = kurtosis_x_weight,
    'min_x_weight' = min_x_weight,
    'max_x_weight' = max_x_weight,
    'range_x_weight' = range_x_weight,
    'Q10_weight' = Q10_weight,
    'Q20_weight' = Q20_weight,
    'Q25_weight' = Q25_weight,
    'Q30_weight' = Q30_weight,
    'Q40_weight' = Q40_weight,
    'Q60_weight' = Q60_weight,
    'Q70_weight' = Q70_weight,
    'Q75_weight' = Q75_weight,
    'Q80_weight' = Q80_weight,
    'Q90_weight' = Q90_weight,
    'theil_index_weight' = theil_index_weight,
    'gini_coeff_weight' = gini_coeff_weight,
    'mean_x_closeness' = mean_x_closeness,
    'median_x_closeness' = median_x_closeness,
    'sd_x_closeness' = sd_x_closeness,
    'iqr_x_closeness' = iqr_x_closeness,
    'skewness_x_closeness' = skewness_x_closeness,
    'kurtosis_x_closeness' = kurtosis_x_closeness,
    'min_x_closeness' = min_x_closeness,
    'max_x_closeness' = max_x_closeness,
    'range_x_closeness' = range_x_closeness,
    'Q10_closeness' = Q10_closeness,
    'Q20_closeness' = Q20_closeness,
    'Q25_closeness' = Q25_closeness,
    'Q30_closeness' = Q30_closeness,
    'Q40_closeness' = Q40_closeness,
    'Q60_closeness' = Q60_closeness,
    'Q70_closeness' = Q70_closeness,
    'Q75_closeness' = Q75_closeness,
    'Q80_closeness' = Q80_closeness,
    'Q90_closeness' = Q90_closeness,
    'theil_index_closeness' = theil_index_closeness,
    'gini_coeff_closeness' = gini_coeff_closeness,
    'mean_x_betweenness' = mean_x_betweenness,
    'median_x_betweenness' = median_x_betweenness,
    'sd_x_betweenness' = sd_x_betweenness,
    'iqr_x_betweenness' = iqr_x_betweenness,
    'skewness_x_betweenness' = skewness_x_betweenness,
    'kurtosis_x_betweenness' = kurtosis_x_betweenness,
    'min_x_betweenness' = min_x_betweenness,
    'max_x_betweenness' = max_x_betweenness,
    'range_x_betweenness' = range_x_betweenness,
    'Q10_betweenness' = Q10_betweenness,
    'Q20_betweenness' = Q20_betweenness,
    'Q25_betweenness' = Q25_betweenness,
    'Q30_betweenness' = Q30_betweenness,
    'Q40_betweenness' = Q40_betweenness,
    'Q60_betweenness' = Q60_betweenness,
    'Q70_betweenness' = Q70_betweenness,
    'Q75_betweenness' = Q75_betweenness,
    'Q80_betweenness' = Q80_betweenness,
    'Q90_betweenness' = Q90_betweenness,
    'theil_index_betweenness' = theil_index_betweenness,
    'gini_coeff_betweenness' = gini_coeff_betweenness,
    'mean_x_degree' = mean_x_degree,
    'median_x_degree' = median_x_degree,
    'sd_x_degree' = sd_x_degree,
    'iqr_x_degree' = iqr_x_degree,
    'skewness_x_degree' = skewness_x_degree,
    'kurtosis_x_degree' = kurtosis_x_degree,
    'min_x_degree' = min_x_degree,
    'max_x_degree' = max_x_degree,
    'range_x_degree' = range_x_degree,
    'Q10_degree' = Q10_degree,
    'Q20_degree' = Q20_degree,
    'Q25_degree' = Q25_degree,
    'Q30_degree' = Q30_degree,
    'Q40_degree' = Q40_degree,
    'Q60_degree' = Q60_degree,
    'Q70_degree' = Q70_degree,
    'Q75_degree' = Q75_degree,
    'Q80_degree' = Q80_degree,
    'Q90_degree' = Q90_degree,
    'theil_index_degree' = theil_index_degree,
    'gini_coeff_degree' = gini_coeff_degree,
    'mean_x_triangle_area' = mean_x_triangle_area,
    'median_x_triangle_area' = median_x_triangle_area,
    'sd_x_triangle_area' = sd_x_triangle_area,
    'iqr_x_triangle_area' = iqr_x_triangle_area,
    'skewness_x_triangle_area' = skewness_x_triangle_area,
    'kurtosis_x_triangle_area' = kurtosis_x_triangle_area,
    'min_x_triangle_area' = min_x_triangle_area,
    'max_x_triangle_area' = max_x_triangle_area,
    'range_x_triangle_area' = range_x_triangle_area,
    'Q10_triangle_area' = Q10_triangle_area,
    'Q20_triangle_area' = Q20_triangle_area,
    'Q25_triangle_area' = Q25_triangle_area,
    'Q30_triangle_area' = Q30_triangle_area,
    'Q40_triangle_area' = Q40_triangle_area,
    'Q60_triangle_area' = Q60_triangle_area,
    'Q70_triangle_area' = Q70_triangle_area,
    'Q75_triangle_area' = Q75_triangle_area,
    'Q80_triangle_area' = Q80_triangle_area,
    'Q90_triangle_area' = Q90_triangle_area,
    'theil_index_triangle_area' = theil_index_triangle_area,
    'gini_coeff_triangle_area' = gini_coeff_triangle_area,
    'mean_x_triangle_perimeter' = mean_x_triangle_perimeter,
    'median_x_triangle_perimeter' = median_x_triangle_perimeter,
    'sd_x_triangle_perimeter' = sd_x_triangle_perimeter,
    'iqr_x_triangle_perimeter' = iqr_x_triangle_perimeter,
    'skewness_x_triangle_perimeter' = skewness_x_triangle_perimeter,
    'kurtosis_x_triangle_perimeter' = kurtosis_x_triangle_perimeter,
    'min_x_triangle_perimeter' = min_x_triangle_perimeter,
    'max_x_triangle_perimeter' = max_x_triangle_perimeter,
    'range_x_triangle_perimeter' = range_x_triangle_perimeter,
    'Q10_triangle_perimeter' = Q10_triangle_perimeter,
    'Q20_triangle_perimeter' = Q20_triangle_perimeter,
    'Q25_triangle_perimeter' = Q25_triangle_perimeter,
    'Q30_triangle_perimeter' = Q30_triangle_perimeter,
    'Q40_triangle_perimeter' = Q40_triangle_perimeter,
    'Q60_triangle_perimeter' = Q60_triangle_perimeter,
    'Q70_triangle_perimeter' = Q70_triangle_perimeter,
    'Q75_triangle_perimeter' = Q75_triangle_perimeter,
    'Q80_triangle_perimeter' = Q80_triangle_perimeter,
    'Q90_triangle_perimeter' = Q90_triangle_perimeter,
    'theil_index_triangle_perimeter' = theil_index_triangle_perimeter,
    'gini_coeff_triangle_perimeter' = gini_coeff_triangle_perimeter
  ))


}
