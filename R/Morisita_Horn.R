#' Morisita_Horn index
#'
#' This is the function that calculates the Morisita_Horn index
#'
#' @param x The coordinate x
#' @param y The coordinate y
#'
#' @return This function returns the Morisita_Horn index
#'
#' @author Xiao Li, \email{xiao.li.xl2@roche.com}
#'
#' @export




Morisita_Horn <-  function(x, y){

  num = 2*sum(x*y)
  den = sum(x^2) + sum(y^2)
  mh = num/den

  return(mh)
}
