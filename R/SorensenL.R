#' Sorensen's L index
#'
#' This is the function that calculates the Sorensen's L index
#'
#' @param x The coordinate x
#' @param y The coordinate y
#'
#' @return This function returns the Sorensen's L index
#'
#' @author Xiao Li, \email{li.xiao@gene.com}
#'
#' @export




SorensenL <- function(x, y){
  num = 0
  for (i in 1:length(x)) {
    num = num + min(x[i], y[i])
  }

  den = sum(x) + sum(y)
  SL = 2*num/den

  return(SL)
}
