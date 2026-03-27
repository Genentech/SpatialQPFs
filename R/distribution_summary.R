#' distribution_summary
#'
#' This is the function that calculates the summary features from a distribution
#'
#' @param x The input vector
#'
#' @return This function returns the summary features from a distribution
#'
#'
#' @importFrom moments skewness kurtosis
#' @importFrom ineq Theil Gini
#'
#' @author Xiao Li, \email{xiao.li.xl2@roche.com}
#'
#' @export


distribution_summary <- function(x){
  
  x = x[!is.na(x)]
  
  # ---- 1. Central Tendency ----
  mean_x <- mean(x)
  median_x <- median(x)
  
  
  # ---- 2. Spread ----
  sd_x <- sd(x)
  iqr_x <- IQR(x)
  
  # ---- 3. Shape ----
  skewness_x <- moments::skewness(x) # Requires `moments` package
  kurtosis_x <- moments::kurtosis(x)
  
  # ---- 4. Extremes ----
  min_x <- min(x)
  max_x <- max(x)
  range_x <- max_x - min_x
  
  
  # ---- 5. Quantiles ----
  quantiles <- quantile(x, probs = c(0.1, 0.2, 0.25, 0.3, 0.4, 0.6, 0.7, 0.75, 0.8, 0.9))
  
  quantile_df <- data.frame(t(quantiles))
  colnames(quantile_df) <- c("Q10", "Q20", "Q25", "Q30", "Q40", "Q60", "Q70", "Q75", "Q80", "Q90")
  
  # ---- 6. Heterogeneity ----
  # Theil Index
  normalized_x <- (x - min(x)) / (max(x) - min(x))
  theil_index <- ineq::Theil(normalized_x)
  
  # Gini Coefficient
  gini_coeff <- ineq::Gini(x)
  
  
  return(data.frame(mean_x, median_x, 
                    sd_x, iqr_x,
                    skewness_x, kurtosis_x,
                    min_x, max_x, range_x,
                    quantile_df,
                    theil_index,
                    gini_coeff))
  
  
}




