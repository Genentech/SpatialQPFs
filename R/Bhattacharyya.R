#' A Bhattacharyya function
#'
#' This is the function that calculates the Bhattacharyya's coefficient index
#'
#' @param x The coordinate x
#' @param y The coordinate y
#'
#' @return This function returns the Bhattacharyya's coefficient index
#'
#' @author Xiao Li, \email{li.xiao@gene.com}
#'
#' @export



Bhattacharyya <- function(x, y){

  BC = sum(sqrt(x*y))

  return(BC)
}
