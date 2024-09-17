#' Bivariate Moran's I index
#'
#' This is the function that calculates the Bivariate Moran's I feature
#'
#' @param x The coordinate x
#' @param y The coordinate y
#' @param W The spatial weight matrix
#'
#' @return This function returns the Bivariate Moran's I feature
#'
#' @author Xiao Li, \email{xiao.li.xl2@roche.com}
#'
#' @export




moran_I_BV <- function(x, y, W){
  # if(is.null(y)) y = x

  xp <- scale(x)[, 1]
  yp <- scale(y)[, 1]
  W[which(is.na(W))] <- 0
  n <- nrow(W)

  global <- (n/sum(W))*(xp%*%W%*%yp)/(n - 1)
  # ## W should be row-standardized
  # W_s <- t(apply(W, 1, function(x) {
  #                         if (sum(x) > 0) {
  #                           x/sum(x)
  #                         } else {
  #                           x
  #                         }
  #                       }
  #         )
  #       )
  # global <- (xp%*%W_s%*%yp)/(n - 1)
  # # local  <- (xp*W_s%*%yp)

  #list(global = global, local  = as.numeric(local))
  return("global" = global)
}
