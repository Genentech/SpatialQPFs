#' A regionalization_crosscorrelation_lattice function
#'
#' This is the core function that calculates the features for Lee's L statistic, using spatial lattice process.
#'
#' @param spp_df The input spatial data.frame, need to have 3 columns: "x_coord", "y_coord" and "cell_class"
#' @param from_type The cell type one wants to use as the "from" cell type
#' @param to_type The cell type one wants to use as the "to" cell type
#' @param scale The spatial range that user wants to investigate
#' @param myplot Whether to plot the results, if available, by default it is set as FALSE
#'
#' @return This core function calculates the features for Lee's L statistic, both global and local 
#'
#' @author Xiao Li, \email{li.xiao@gene.com}
#'
#' @export





regionalization_crosscorrelation_lattice <- function(spp_df, from_type, to_type, scale, myplot){

  L = max(max(spp_df$x_coord), max(spp_df$y_coord))
  radius = scale/L


  spp_tumor = spp_df %>% filter(cell_class == to_type) %>% select(x, y)
  spp_immune = spp_df %>% filter(cell_class == from_type) %>% select(x, y)

  # r <- raster(ncol=50, nrow=50, xmn = 0, xmx = 1, ymn = 0, ymx = 1)
  xg <- seq(0, max(spp_df$x), length.out = ceiling(max(spp_df$x)/(2*radius))+1)
  yg <- seq(0, max(spp_df$y), length.out = ceiling(max(spp_df$y)/(2*radius))+1)

  grid <- expand.grid(xg, yg)


  binxy_TC <- data.frame(x=findInterval(spp_tumor$x, xg, rightmost.closed = TRUE),
                         y=findInterval(spp_tumor$y, yg, rightmost.closed = TRUE))

  binxy_IC <- data.frame(x=findInterval(spp_immune$x, xg, rightmost.closed = TRUE),
                         y=findInterval(spp_immune$y, yg, rightmost.closed = TRUE))

  results_TC <- table(binxy_TC)
  results_IC <- table(binxy_IC)


  d_TC <- as.data.frame.table(results_TC)
  d_IC <- as.data.frame.table(results_IC)

  d_TC$x = as.numeric(as.character(d_TC$x))
  d_IC$x = as.numeric(as.character(d_IC$x))
  d_TC$y = as.numeric(as.character(d_TC$y))
  d_IC$y = as.numeric(as.character(d_IC$y))



  xx <- xg[-length(xg)] + 0.5*diff(xg)
  d_TC$xg <- xx[d_TC$x]
  d_IC$xg <- xx[d_IC$x]
  yy <- yg[-length(yg)] + 0.5*diff(yg)
  d_TC$yg <- yy[d_TC$y]
  d_IC$yg <- yy[d_IC$y]


  ### make a background grid space

  bk_grid = cbind(expand.grid(1:length(xx), 1:length(yy)), expand.grid(xx, yy))
  names(bk_grid) = c('x', 'y', 'xg', 'yg')


  d_TC = suppressMessages(left_join(bk_grid, d_TC))
  d_IC = suppressMessages(left_join(bk_grid, d_IC))

  
  ROI = !(is.na(d_TC$Freq) | is.na(d_IC$Freq) | (d_TC$Freq == 0 & d_IC$Freq == 0))
  
  d_TC$p = d_TC$Freq/nrow(spp_tumor)
  d_IC$p = d_IC$Freq/nrow(spp_immune)
  
  
  d_TC_poly = d_TC %>% filter(ROI)
  d_IC_poly = d_IC %>% filter(ROI)
  
  
  d_comb_poly = suppressMessages(left_join(d_TC_poly %>% rename(p_TC = p), 
                                           d_IC_poly %>% rename(p_IC = p),
                                           by = c("xg", "yg")) )
  
  coordinates(d_comb_poly) = ~ xg + yg
  
  
  nb_comb <- dnearneigh(d_comb_poly, d1 = 0, d2 = 2*radius + 0.0001)
  
  
  W_comb <- nb2listw(nb_comb, style="B", zero.policy=TRUE)
  
  simula_lee <- function(x, y, listw, nsim = nsim, zero.policy = NULL, na.action = na.fail) {
    
    if (deparse(substitute(na.action)) == "na.pass") 
      stop ("na.pass not permitted")
    na.act <- attr(na.action(cbind(x, y)), "na.action")
    x[na.act] <- NA
    y[na.act] <- NA
    x <- na.action(x)
    y <- na.action(y)
    if (!is.null(na.act)) {
      subset <- !(1:length(listw$neighbours) %in% na.act)
      listw <- subset(listw, subset, zero.policy = zero.policy)
    }
    n <- length(listw$neighbours)
    if ((n != length(x)) | (n != length(y))) 
      stop ("objects of different length")
    gamres <- suppressWarnings(nsim > gamma(n + 1))
    if (gamres) 
      stop ("nsim too large for this number of observations")
    if (nsim < 1) 
      stop ("nsim too small")
    xy <- data.frame(x, y)
    S2 <- sum((unlist(lapply(listw$weights, sum)))^2)
    
    lee_boot <- function(var, i, ...) {
      return(spdep::lee(x = var[i, 1], y = var[i, 2], ...)$localL)
    }
    
    res <- boot::boot(xy, statistic = lee_boot, R = nsim, sim = "permutation", 
                      listw = listw, n = n, S2 = S2, zero.policy = zero.policy)
  }
  
  
  
  ########################################################################################################
  ############################        from from_type to to_type           ################################
  ########################################################################################################
  
  Lee_IC_TC <- spdep::lee(d_comb_poly$p_IC, d_comb_poly$p_TC, W_comb, zero.policy=TRUE, length(nb_comb), Szero(W_comb), NAOK=TRUE)
  
  
  local_sims_IC_TC <- simula_lee(d_comb_poly$p_IC,
                                 d_comb_poly$p_TC,
                                 W_comb,
                                 nsim = 100,
                                 zero.policy = TRUE,
                                 na.action = na.omit)
  
  
  m_i <- Lee_IC_TC$localL  # local values
  
  # Identify the significant values 
  alpha <- 0.05  # for a 95% confidence interval
  probs <- c(alpha/2, 1-alpha/2)
  intervals <- t(apply(t(local_sims_IC_TC[[2]]), 1, function(x) quantile(x, probs = probs)))
  sig <- ( m_i < intervals[ , 1] ) | ( m_i > intervals[ , 2] )
  
  
  
  # Prepare for plotting
  map_sf_IC_TC <- sf::st_as_sf(data.frame(d_comb_poly), coords = c("xg", "yg"))
  map_sf_IC_TC$sig <- sig
  
  # Identify the Lee's L clusters
  Xp <- scale(d_comb_poly$p_IC)[,1]
  Yp <- scale(d_comb_poly$p_TC)[,1]
  
  # Adjacency Matrix
  
  W  <- nb2mat(nb_comb, style="W", zero.policy=TRUE)
  W[which(is.na(W))] <- 0
  
  
  patterns <- as.character(interaction(Xp > 0, W %*% Yp > 0)) 
  patterns <- patterns %>% 
    stringr::str_replace_all("TRUE","High") %>% 
    stringr::str_replace_all("FALSE","Low")
  patterns[map_sf_IC_TC$sig == 0] <- "non-significant"
  map_sf_IC_TC$patterns <- patterns
  
  # Rename Lee's L clusters
  map_sf_IC_TC$local_pattern <- factor(map_sf_IC_TC$patterns,
                                       levels = c("High.High", "High.Low", "Low.High", "Low.Low", "non-significant"),
                                       labels = c("High-High", "High-Low", "Low-High", "Low-Low", "non-significant"))
  
  
  
  if (myplot){
    print(ggplot() + 
            geom_sf(data = map_sf_IC_TC, aes(color = local_pattern)) +
            labs(color = "Local Lee's L clusters") + xlab("") + ylab("") +
            scale_colour_manual(values = c("red", "blue", "green", "orange", "white")) +
            theme(axis.text.x=element_blank(),
                  axis.ticks.x=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks.y=element_blank()
            )  +
            ggtitle(paste0("Local Lee's L from ", from_type, " to ", to_type)) + 
            theme(plot.title = element_text(hjust = 0.5)))
  }
  
  
  
  
  Lee_HL_IC_TC = sum(map_sf_IC_TC$local_pattern == "High-Low")/length(map_sf_IC_TC$local_pattern)
  Lee_HH_IC_TC = sum(map_sf_IC_TC$local_pattern == "High-High")/length(map_sf_IC_TC$local_pattern)
  Lee_LH_IC_TC = sum(map_sf_IC_TC$local_pattern == "Low-High")/length(map_sf_IC_TC$local_pattern)
  Lee_LL_IC_TC = sum(map_sf_IC_TC$local_pattern == "Low-Low")/length(map_sf_IC_TC$local_pattern)
  
  

  
  ########################################################################################################
  ############################        from to_type to from_type           ################################
  ########################################################################################################
  
  Lee_TC_IC <- spdep::lee(d_comb_poly$p_TC, d_comb_poly$p_IC, W_comb, zero.policy=TRUE, length(nb_comb), Szero(W_comb), NAOK=TRUE)
  
  
  local_sims_TC_IC <- simula_lee(d_comb_poly$p_TC,
                                 d_comb_poly$p_IC,
                                 W_comb,
                                 nsim = 100,
                                 zero.policy = TRUE,
                                 na.action = na.omit)
  
  
  m_i <- Lee_TC_IC$localL  # local values
  
  # Identify the significant values 
  alpha <- 0.05  # for a 95% confidence interval
  probs <- c(alpha/2, 1-alpha/2)
  intervals <- t(apply(t(local_sims_TC_IC[[2]]), 1, function(x) quantile(x, probs = probs)))
  sig <- ( m_i < intervals[ , 1] ) | ( m_i > intervals[ , 2] )
  
  
  
  # Prepare for plotting
  map_sf_TC_IC <- sf::st_as_sf(data.frame(d_comb_poly), coords = c("xg", "yg"))
  map_sf_TC_IC$sig <- sig
  
  # Identify the Lee's L clusters
  Xp <- scale(d_comb_poly$p_TC)[,1]
  Yp <- scale(d_comb_poly$p_IC)[,1]
  
  # Adjacency Matrix
  
  W  <- nb2mat(nb_comb, style="W", zero.policy=TRUE)
  W[which(is.na(W))] <- 0
  
  
  patterns <- as.character(interaction(Xp > 0, W %*% Yp > 0)) 
  patterns <- patterns %>% 
    stringr::str_replace_all("TRUE","High") %>% 
    stringr::str_replace_all("FALSE","Low")
  patterns[map_sf_TC_IC$sig == 0] <- "non-significant"
  map_sf_TC_IC$patterns <- patterns
  
  # Rename Lee's L clusters
  map_sf_TC_IC$local_pattern <- factor(map_sf_TC_IC$patterns,
                                       levels = c("High.High", "High.Low", "Low.High", "Low.Low", "non-significant"),
                                       labels = c("High-High", "High-Low", "Low-High", "Low-Low", "non-significant"))
  
  
  
  if (myplot){
    print(ggplot() + 
            geom_sf(data = map_sf_TC_IC, aes(color = local_pattern)) +
            labs(color = "Local Lee's L clusters") + xlab("") + ylab("") +
            scale_colour_manual(values = c("red", "blue", "green", "orange", "white")) +
            theme(axis.text.x=element_blank(),
                  axis.ticks.x=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks.y=element_blank()
            )  +
            ggtitle(paste0("Local Lee's L from ", to_type, " to ", from_type)) + 
            theme(plot.title = element_text(hjust = 0.5)))
  }
  
  
  
  
  Lee_HL_TC_IC = sum(map_sf_TC_IC$local_pattern == "High-Low")/length(map_sf_TC_IC$local_pattern)
  Lee_HH_TC_IC = sum(map_sf_TC_IC$local_pattern == "High-High")/length(map_sf_TC_IC$local_pattern)
  Lee_LH_TC_IC = sum(map_sf_TC_IC$local_pattern == "Low-High")/length(map_sf_TC_IC$local_pattern)
  Lee_LL_TC_IC = sum(map_sf_TC_IC$local_pattern == "Low-Low")/length(map_sf_TC_IC$local_pattern)
  
  
  

  return(list("Lee_L" = Lee_TC_IC$L,
              "Lee_HL_TC_IC" = Lee_HL_TC_IC, "Lee_HH_TC_IC" = Lee_HH_TC_IC, "Lee_LH_TC_IC" = Lee_LH_TC_IC, "Lee_LL_TC_IC" = Lee_LL_TC_IC,
              "Lee_HL_IC_TC" = Lee_HL_IC_TC, "Lee_HH_IC_TC" = Lee_HH_IC_TC, "Lee_LH_IC_TC" = Lee_LH_IC_TC, "Lee_LL_IC_TC" = Lee_LL_IC_TC
             
  ))
}

