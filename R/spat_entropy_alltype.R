#' Main function to generate spatial entropy features for all cell types involves
#'
#' This is a function that calculates the spatial entropy when all cell types are accounted. The methods included are Shannon entropy, co-occurrence based Shannon entropy, Altieri entropy, Leibovici entropy and Claramunt entropy (spatial form of diversity)
#'
#' @param path The path for the directory that contains the data
#' @param file The file name in the directory
#' @param cell_class_var The column name in each file that indicates the cell class variable
#' @param x_var The column name in each file that indicates the cell location x variable
#' @param y_var The column name in each file that indicates the cell location y variable
#' @param scale The spatial range that user wants to investigate
#' @param side_length The side length of squared FOVs
#' @param num_FOV Number of FOVs to choose
#' @param set_seed set seed for reproducibility
#' @param myplot Whether to plot the results, if available, by default it is set as FALSE
#'
#' @return This function returns different version of spatial entropy
#' \item{shannon}{shannon entropy}
#' \item{shannonZ_entropy}{co-occurrence based Shannon entropy}
#' \item{altieri_entropy_SMI}{Spatial mutual information of Altieri entropy}
#' \item{altieri_entropy_RES}{Global spatial residual of Altieri entropy}
#' \item{leibovici_entropy}{Leibovici entropy at specified distances}
#' \item{spat_diversity}{Spatial form of diversity, i.e. Claramunt entropy}
#'
#' @import SpatEntropy
#' @import dplyr
#'
#' @author Xiao Li, \email{xiao.li.xl2@roche.com}
#'
#' @export




spat_entropy_alltype <- function(path = '/Users/lix233/Haystack/5813_cell_centers/',
                                 file = '267a5ed2-b6ee-4973-bf2d-3526b668eda0.csv',
                                 cell_class_var = "cell_class",
                                 x_var = "X0", y_var = "X1",
                                 scale = 200,
                                 side_length = 2000, num_FOV = 5, set_seed = 42, myplot = FALSE){

  spp_df = read.csv(paste0(path, file)) %>%
    select(x_coord = all_of(x_var), y_coord = all_of(y_var),
           cell_class = all_of(cell_class_var)) %>%
    mutate(y_coord = max(y_coord)-y_coord)
  
  
  if (((max(spp_df$x) - min(spp_df$x)) <= side_length) | ((max(spp_df$y) - min(spp_df$y)) <= side_length)) {
    
    L = max(max(spp_df$x_coord)-min(spp_df$x_coord), max(spp_df$y_coord)-min(spp_df$y_coord))
    
    
    spp_df = spp_df %>% mutate(x = normalized_coords(spp_df$x_coord - min(spp_df$x_coord), L),
                               y = normalized_coords(spp_df$y_coord - min(spp_df$y_coord), L),
                               cell_class = as.factor(spp_df$cell_class))
    
    shannon_entropy = SpatEntropy::shannon(as.vector(spp_df$cell_class))
    
    W = owin(xrange = c(min(spp_df$x), max(spp_df$x)),
             yrange = c(min(spp_df$y), max(spp_df$y)))
    dat_ppp = ppp(x = spp_df$x, y = spp_df$y, marks = spp_df$cell_class, window = W)
    
    radius = scale/L
    
    altieri_entropy = altieri(data = dat_ppp, distbreak = c(radius, 2*radius), verbose = F, plotout = myplot)
    
    leibovici_entropy = leibovici(dat_ppp, ccdist = radius, plotout = myplot)
    
    spat_diversity = spatial_diversity(spp_df, W)
    
    res_list = list(
      "shannon" = shannon_entropy$shann,
      "shannonZ_entropy" = altieri_entropy$ShannonZ$shannZ,
      "altieri_entropy_SMI" = altieri_entropy$SMI,
      "altieri_entropy_RES" = altieri_entropy$RES,
      "leibovici_entropy" = leibovici_entropy$leib,
      "spat_diversity" = spat_diversity
    )
    
    
  } else {
    ## random crop the matrix (num_FOV) times and median the feature vector
    
    D_res = data.frame();
    spp_df = spp_df %>% rename(x = x_coord, y = y_coord)
    
    for (repeat_FOV in 1:num_FOV){
      print(repeat_FOV)
      x_origin = round(runif(1, min = min(spp_df$x), max = max(spp_df$x)-side_length))
      y_origin = round(runif(1, min = min(spp_df$y), max = max(spp_df$y)-side_length))
      
      spp_df_sub = spp_df %>% 
        filter((x >= x_origin) & (x < (x_origin + side_length)) & (y >= y_origin) & (y < (y_origin + side_length)))
      
      c = 0
      
      while(length(unique(spp_df_sub$cell_class)) <= 1 & nrow(spp_df_sub) <= 10 & c < 1000) {
        
        x_origin = round(runif(1, min = min(spp_df$x), max = max(spp_df$x)-side_length))
        y_origin = round(runif(1, min = min(spp_df$y), max = max(spp_df$y)-side_length))
        
        spp_df_sub = spp_df %>%
          filter((x >= x_origin) & (x < (x_origin + side_length)) & (y >= y_origin) & (y < (y_origin + side_length)))
        # count repeated times
        c = c+1
      }
      
      if (c == 1000){
        next
      } else {
        
        shannon_entropy = SpatEntropy::shannon(as.vector(spp_df_sub$cell_class))
        
        W = owin(xrange = c(min(spp_df_sub$x), max(spp_df_sub$x)),
                 yrange = c(min(spp_df_sub$y), max(spp_df_sub$y)))
        dat_ppp = ppp(x = spp_df_sub$x, y = spp_df_sub$y, marks = spp_df_sub$cell_class, window = W)
        
        altieri_entropy = altieri(data = dat_ppp, distbreak = c(scale, 2*scale), verbose = F, plotout = myplot)
        
        leibovici_entropy = leibovici(dat_ppp, ccdist = scale, plotout = myplot)
        
        spat_diversity = spatial_diversity(spp_df, W)
        
        
        this_res = data.frame("shannon" = shannon_entropy$shann,
                              "shannonZ_entropy" = altieri_entropy$ShannonZ$shannZ,
                              "altieri_entropy_SMI" = altieri_entropy$SMI,
                              "altieri_entropy_RES" = altieri_entropy$RES,
                              "leibovici_entropy" = leibovici_entropy$leib,
                              "spat_diversity" = spat_diversity)
        D_res = rbind(D_res, this_res)
        
      }
      
    }
    
    if (dim(D_res)[1] > 0) {
      D_res = data.frame(t(apply(D_res, 2, function(x) median(x[which(!is.na(x))]))))
      
      
      res_list = list(
        "shannon" = D_res$shannon,
        "shannonZ_entropy" = D_res$shannonZ_entropy,
        "altieri_entropy_SMI" = D_res$altieri_entropy_SMI,
        "altieri_entropy_RES" = D_res$altieri_entropy_RES,
        "leibovici_entropy" = D_res$leibovici_entropy,
        "spat_diversity" = D_res$spat_diversity
      )
    } else {
      res_list = list(
        "shannon" = NA,
        "shannonZ_entropy" = NA,
        "altieri_entropy_SMI" = NA,
        "altieri_entropy_RES" = NA,
        "leibovici_entropy" = NA,
        "spat_diversity" = NA
      )
    }
    
  }
  


  return(res_list)

}
