#' Main function to generate local entropy for each single cell 
#' 
#' This is a function that calculates the local entropy for each single cell 
#' 
#' @param path The path for the directory that contains the data
#' @param file The file name in the directory
#' @param cell_class_var The column name in each file that indicates the cell class variable
#' @param x_var The column name in each file that indicates the cell location x variable
#' @param y_var The column name in each file that indicates the cell location y variable
#' @param scale The spatial range that user wants to define as local neighborhood
#' @param myplot Whether to plot the results, if available, by default it is set as FALSE
#'
#' @return This function returns a vector of local entropy for each single cell 
#' \item{local_entropy}{a vector of local entropy for each single cell}
#'
#' @importFrom RANN nn2
#' @import dplyr
#'
#' @author Xiao Li, \email{xiao.li.xl2@roche.com}
#'
#' @export






local_entropy <- function(path = '/Users/lix233/Haystack/5813_cell_centers/',
                          file = '267a5ed2-b6ee-4973-bf2d-3526b668eda0.csv',
                          cell_class_var = "cell_class",
                          x_var = "X0", y_var = "X1",
                          scale = 200,  
                          myplot = T){
  
  spp_df = read.csv(paste0(path, file)) %>%
    select(x_coord = all_of(x_var), y_coord = all_of(y_var),
           cell_class = all_of(cell_class_var)) %>%
    mutate(y_coord = max(y_coord)-y_coord)
  
  
  target_spp_df = spp_df
  
  local_entropy = c()
  for (i in 1:nrow(target_spp_df)){
    
    this_point = target_spp_df[i, ]
    
    this_nb = RANN::nn2(data = this_point[, c('x_coord','y_coord')], query = spp_df[, c('x_coord','y_coord')], searchtype = "radius", radius = scale)
    
    
    idx = which(this_nb$nn.idx > 0) ## all the nb idx, including itself
    idx = idx[idx != which(this_nb$nn.dists == 0)] ## exclude itself: the distance between the point to itself is exactly 0 !
    
    if (length(idx) > 0){
      p = table(spp_df[idx, "cell_class"])/length(idx)
      
      local_entropy = c(local_entropy, sum(-p*log2(p)))
      
    } else {
      local_entropy = c(local_entropy, NA)
    }
    
  }
  
  if(myplot){
    
    this_plot = ggplot(target_spp_df) + aes (x = x_coord, y = y_coord, colour = local_entropy) +
      geom_point() +
      scale_colour_gradientn(colours = rainbow(n = 10, rev = T)) 
    print(this_plot)
    
  }
  
  
  return("local_entropy" = local_entropy)
  
  
}
