#' Main function to generate spatial entropy features for one cell type
#'
#' This is a function that calculates the spatial entropy when all cell types are accounted. The methods included are batty entropy and battyLISA entropy (Karlstrom and Ceccato's entropy)
#'
#' @param path The path for the directory that contains the data
#' @param file The file name in the directory
#' @param cell_class_var The column name in each file that indicates the cell class variable
#' @param x_var The column name in each file that indicates the cell location x variable
#' @param y_var The column name in each file that indicates the cell location y variable
#' @param cell_type The cell type one wants to use for calculating spatial entropy features
#' @param scale The sides of the square grids
#' @param myplot Whether to plot the results, if available, by default it is set as FALSE
#'
#' @return This function returns different version of spatial entropy
#' \item{batty_entropy}{batty entropy}
#' \item{battyLISA_entropy}{batty's LISA entropy, i.e. Karlstrom and Ceccato's entropy}
#'
#' @import SpatEntropy
#' @import dplyr
#'
#' @author Xiao Li, \email{xiao.li.xl2@roche.com}
#'
#' @export




spat_entropy_onetype <- function(path = '/Users/lix233/Haystack/5813_cell_centers/',
                                 file = '267a5ed2-b6ee-4973-bf2d-3526b668eda0.csv',
                                 cell_class_var = "cell_class",
                                 x_var = "X0", y_var = "X1",
                                 cell_type = "Tumor",
                                 scale = 200, myplot = FALSE){

  spp_df = read.csv(paste0(path, file)) %>% filter(get(cell_class_var) == cell_type) %>%
    select(x_coord = all_of(x_var), y_coord = all_of(y_var),
           cell_class = all_of(cell_class_var)) %>%
    mutate(y_coord = max(y_coord)-y_coord)

  L = max(max(spp_df$x_coord)-min(spp_df$x_coord), max(spp_df$y_coord)-min(spp_df$y_coord))

  spp_df = spp_df %>% mutate(x = normalized_coords(spp_df$x_coord - min(spp_df$x_coord), L),
                             y = normalized_coords(spp_df$y_coord - min(spp_df$y_coord), L),
                             cell_class = as.factor(spp_df$cell_class))

  W = owin(xrange = c(min(spp_df$x), max(spp_df$x)),
           yrange = c(min(spp_df$y), max(spp_df$y)))
  dat_ppp = ppp(x = spp_df$x, y = spp_df$y, marks = spp_df$cell_class, window = W)

  radius = scale/L

  xg <- radius + 2*radius*(0:(floor(max(spp_df$x)/(2*radius))-1))
  yg <- radius + 2*radius*(0:(floor(max(spp_df$y)/(2*radius))-1))

  cen = expand.grid(xg, yg)

  # total_area = area(W)
  # num_part = round(total_area/(2*radius)^2)

  batty_entropy = batty(dat_ppp, partition = as.matrix(cen), category = cell_type, plotout = myplot)
  battyLISA_entropy = battyLISA(dat_ppp, partition = as.matrix(cen), category = cell_type, plotout = myplot)


  return(list("batty_entropy" = batty_entropy$batty,
              "battyLISA_entropy" = battyLISA_entropy$karlstrom))

}
