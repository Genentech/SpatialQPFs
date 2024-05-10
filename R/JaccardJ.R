#' A JaccardJ function
#'
#' This is the function that calculates the Jaccard's J index
#'
#' @param x The coordinate x
#' @param y The coordinate y
#'
#' @return This function returns the Jaccard's J index
#'
#' @author Xiao Li, \email{li.xiao@gene.com}
#'
#' @export


JaccardJ <- function(x, y){
  num = 0
  for (i in 1:length(x)) {
    num = num + min(x[i], y[i])
  }

  den = sum(x) + sum(y) - num
  JJ = num/den

  return(JJ)
}
