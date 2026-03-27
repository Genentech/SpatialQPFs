#' Spatial Form of Diversity
#'
#' This is a function that computes the spatial version of shannon entropy: https://link.springer.com/chapter/10.1007/11556114_14
#'
#' @param spp_df The data frame that contains the coordinates and cell class columns
#' @param W The study area of the spatial data
#'
#' @return This function returns the spatial form of diversity
#'
#' @author Xiao Li, \email{xiao.li.xl2@roche.com}
#'
#' @export




spatial_diversity <- function(spp_df, W){

  tab = table(spp_df$cell_class)
  p_i = tab/sum(tab)
  total_num_types = length(tab)

  int = ext = c()
  for(i in names(tab)){

    this_class = spp_df %>% filter(cell_class == i) %>% select(x, y)
    other_class = spp_df %>% filter(cell_class != i) %>% select(x, y)

    this_class_ppp = ppp(x = this_class$x, y = this_class$y, window = W)
    other_class_ppp = ppp(x = other_class$x, y = other_class$y, window = W)

    int = c(int, mean(pairdist(this_class_ppp)))
    ext = c(ext, mean(crossdist(this_class_ppp, other_class_ppp)))
  }

  H = -sum(int/ext*p_i*log2(p_i))

  return(H)
}
