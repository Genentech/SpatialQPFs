#' Normalizes the coordinates to 0-1 scale
#'
#' This is a simple function that normalizes the coordinates to 0-1 scale
#'
#' @param x The original coordinate x or y, before normalization
#' @param L The scaling factor for normalization
#'
#' @return This function returns the normalized coordinates
#'
#' @author Xiao Li, \email{li.xiao@gene.com}
#'
#' @export



normalized_coords <- function(x, L) {
  x_prime = x/L
  return(x_prime)
}
