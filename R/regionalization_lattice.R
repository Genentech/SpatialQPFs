#' A regionalization_lattice function
#'
#' This is the core function that calculates the features for areal data, using spatial lattice process.
#'
#' @param spp_df The input spatial data.frame, need to have 3 columns: "x_coord", "y_coord" and "cell_class"
#' @param from_type The cell type one wants to use as the "from" cell type
#' @param to_type The cell type one wants to use as the "to" cell type
#' @param scale The spatial range that user wants to investigate
#' @param myplot Whether to plot the results, if available, by default it is set as FALSE
#'
#' @return This core function returns the features for areal data, using spatial lattice process
#'
#' @author Xiao Li, \email{li.xiao@gene.com}
#'
#' @export





regionalization_lattice <- function(spp_df, from_type, to_type, scale,
                                    myplot){

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

  ### make NA grid to 0
  d_TC$Freq[is.na(d_TC$Freq)] = 0
  d_IC$Freq[is.na(d_IC$Freq)] = 0



  ### make plots
  # plot(grid, t="n", xaxs="i", yaxs="i")
  # points(spp_tumor, col = "black", cex = 0.01)
  # with(d_TC, text(xg, yg, label=Freq))
  #
  # plot(grid, t="n", xaxs="i", yaxs="i")
  # points(spp_immune, col = "grey", cex = 0.01)
  # with(d_IC, text(xg, yg, label=Freq))


  d_TC$p = d_TC$Freq/nrow(spp_tumor)
  d_IC$p = d_IC$Freq/nrow(spp_immune)

  ### Bhattacharyya
  BC = Bhattacharyya(d_TC$p, d_IC$p)

  ### use the spatial lattice method
  mh = Morisita_Horn(d_TC$p, d_IC$p)
  ### Jaccard and Sorensen Indices

  JJ = JaccardJ(d_TC$p, d_IC$p)
  SL = SorensenL(d_TC$p, d_IC$p)


  ### re-assign the d_TC and d_IC
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

  d_TC = suppressMessages(left_join(bk_grid, d_TC))
  d_IC = suppressMessages(left_join(bk_grid, d_IC))


  d_TC$ROI = (!(is.na(d_TC$Freq) | d_TC$Freq == 0)) | (!(is.na(d_IC$Freq) | d_IC$Freq == 0))
  d_IC$ROI = (!(is.na(d_IC$Freq) | d_IC$Freq == 0)) | (!(is.na(d_TC$Freq) | d_TC$Freq == 0))


  d_TC$p = d_TC$Freq/nrow(spp_tumor)
  d_IC$p = d_IC$Freq/nrow(spp_immune)


  d_TC_poly = d_TC %>% filter(ROI)
  d_IC_poly = d_IC %>% filter(ROI)


  ### convert the format into spatialdataframe
  coordinates(d_TC_poly) = ~ xg + yg
  coordinates(d_IC_poly) = ~ xg + yg

  if (any(is.na(d_TC_poly$p))) {
    d_TC_poly$p[is.na(d_TC_poly$p)] = 0
  }

  if (any(is.na(d_IC_poly$p))) {
    d_IC_poly$p[is.na(d_IC_poly$p)] = 0
  }



  ### moran's I for TC and Lymphocytes
  nb_TC <- dnearneigh(d_TC_poly, d1 = 0, d2 = 2*radius + 0.0001)

  W_TC <- nb2listw(nb_TC, style="B", zero.policy=TRUE)

  moran_TC <- spdep::moran(d_TC_poly$p, W_TC, zero.policy=TRUE, length(nb_TC), Szero(W_TC), NAOK=TRUE)$I
  geary_TC <- spdep::geary(x=d_TC_poly$p,listw=W_TC, zero.policy=TRUE, n=length(nb_TC),n1=length(nb_TC)-1,S0= Szero(W_TC))$C

  ### local moran and local geary
  #### local moran's I
  lmoran_TC <- localmoran(d_TC_poly$p, W_TC, zero.policy=TRUE)
  localp = lmoran_TC[, 5]; localp[is.na(localp)] = 1

  moran_HL_TC <- sum(attr(lmoran_TC, 'quadr')$mean == "High-Low" & localp < 0.05)/nrow(attr(lmoran_TC, 'quadr'))
  moran_HH_TC <- sum(attr(lmoran_TC, 'quadr')$mean == "High-High" & localp < 0.05)/nrow(attr(lmoran_TC, 'quadr'))
  moran_LH_TC <- sum(attr(lmoran_TC, 'quadr')$mean == "Low-High" & localp < 0.05)/nrow(attr(lmoran_TC, 'quadr'))
  moran_LL_TC <- sum(attr(lmoran_TC, 'quadr')$mean == "Low-Low" & localp < 0.05)/nrow(attr(lmoran_TC, 'quadr'))

  
  if (myplot){
    localp = lmoran_TC[, 5]; localp[is.na(localp)] = 1

    local_pattern = as.character(attr(lmoran_TC, 'quadr')$mean)
    local_pattern[localp >= 0.05] = "non-significant"

    print(ggplot(d_TC_poly@data, aes(x, y, color = local_pattern)) +
      geom_point(size = 2) +
      labs(color = "Local Moran's I pattern") + xlab("") + ylab("") +
      scale_colour_manual(values = c("red", "blue", "green", "orange", "white")) +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()
      )  +
      ggtitle(paste0("Local Moran's I for ", to_type)) + 
      theme(plot.title = element_text(hjust = 0.5)))
  }
  


  #### local geary's C
  lgeary_TC <- spdep::localC_perm(x = d_TC_poly$p, listw = W_TC, zero.policy=TRUE, nsim = 999)

  localp = attr(lgeary_TC, "pseudo-p")[,5]; localp[is.na(localp)] = 1
  geary_HH_TC <- sum(attr(lgeary_TC, 'cluster') == "High-High" &  localp < 0.05)/length(lgeary_TC)
  geary_LL_TC <- sum(attr(lgeary_TC, 'cluster') == "Low-Low" & localp < 0.05)/length(lgeary_TC)


  if (myplot){
    local_pattern = as.character(attr(lgeary_TC, 'cluster'))
    local_pattern[localp >= 0.05] = "non-significant"

    print(ggplot(d_TC_poly@data, aes(x, y, color = local_pattern)) +
      geom_point(size = 2) +
      labs(color = "Local geary's C pattern") + xlab("") + ylab("") +
      scale_colour_manual(values = c("red", "orange", "green",  "white", "blue")) +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()
      )  +
      ggtitle(paste0("Local Geary's C for ", to_type)) + 
      theme(plot.title = element_text(hjust = 0.5)))
  }


  
  nb_IC <- dnearneigh(d_IC_poly, d1 = 0, d2 = 2*radius + 0.0001)

  W_IC <- nb2listw(nb_IC, style="B", zero.policy=TRUE)

  moran_IC <- spdep::moran(d_IC_poly$p, W_IC, zero.policy=TRUE, length(nb_IC), Szero(W_IC), NAOK=TRUE)$I
  geary_IC <- spdep::geary(x=d_IC_poly$p,listw=W_IC, zero.policy=TRUE, n=length(nb_IC),n1=length(nb_IC)-1,S0= Szero(W_IC))$C

  ### local moran and local geary
  #### local moran's I
  lmoran_IC <- localmoran(d_IC_poly$p, W_IC, zero.policy=TRUE)
  localp = lmoran_IC[, 5]; localp[is.na(localp)] = 1

  moran_HL_IC <- sum(attr(lmoran_IC, 'quadr')$mean == "High-Low" & localp < 0.05)/nrow(attr(lmoran_IC, 'quadr'))
  moran_HH_IC <- sum(attr(lmoran_IC, 'quadr')$mean == "High-High" & localp < 0.05)/nrow(attr(lmoran_IC, 'quadr'))
  moran_LH_IC <- sum(attr(lmoran_IC, 'quadr')$mean == "Low-High" & localp < 0.05)/nrow(attr(lmoran_IC, 'quadr'))
  moran_LL_IC <- sum(attr(lmoran_IC, 'quadr')$mean == "Low-Low" & localp < 0.05)/nrow(attr(lmoran_IC, 'quadr'))

  
  
  if (myplot){
    localp = lmoran_IC[, 5]; localp[is.na(localp)] = 1

    local_pattern = as.character(attr(lmoran_IC, 'quadr')$mean)
    local_pattern[localp >= 0.05] = "non-significant"

    print(ggplot(d_IC_poly@data, aes(x, y, color = local_pattern)) +
      geom_point(size = 2) +
      labs(color = "Local Moran's I pattern") + xlab("") + ylab("") +
      scale_colour_manual(values = c("red", "blue", "green", "orange", "white")) +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()
      )  +
      ggtitle(paste0("Local Moran's I for ", from_type)) + 
      theme(plot.title = element_text(hjust = 0.5)))
  }
  



  #### local geary's C
  lgeary_IC <- spdep::localC_perm(x = d_IC_poly$p, listw = W_IC, zero.policy=TRUE, nsim = 999)

  localp = attr(lgeary_IC, "pseudo-p")[,5]; localp[is.na(localp)] = 1
  geary_HH_IC <- sum(attr(lgeary_IC, 'cluster') == "High-High" &  localp < 0.05)/length(lgeary_IC)
  geary_LL_IC <- sum(attr(lgeary_IC, 'cluster') == "Low-Low" & localp < 0.05)/length(lgeary_IC)


  
  if (myplot){
    local_pattern = as.character(attr(lgeary_IC, 'cluster'))
    local_pattern[localp >= 0.05] = "non-significant"

    print(ggplot(d_IC_poly@data, aes(x, y, color = local_pattern)) +
      geom_point(size = 2) +
      labs(color = "Local geary's C pattern") + xlab("") + ylab("") +
      scale_colour_manual(values = c("red", "orange", "green",  "white", "blue")) +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()
      ) +
      ggtitle(paste0("Local Geary's C for ", from_type)) + 
      theme(plot.title = element_text(hjust = 0.5)))
  }
  





  ### Cross-type (Bivariate) Moran's I

  W  <- nb2mat(nb_TC, style="B", zero.policy=TRUE)

  moran_BV = moran_I_BV(d_TC_poly$p, d_IC_poly$p, W)[1,1]




  ### Getis-Ord methods
  ## include self: to calculate G* which counts itself and more popular

  nb_TC <- dnearneigh(d_TC_poly, d1 = 0, d2 = 2*radius + 0.0001)
  nb_TC <- include.self(nb_TC)
  W_TC <- nb2listw(nb_TC, style="B", zero.policy=TRUE)
  nb_IC <- dnearneigh(d_IC_poly, d1 = 0, d2 = 2*radius + 0.0001)
  nb_IC <- include.self(nb_IC)
  W_IC <- nb2listw(nb_IC, style="B", zero.policy=TRUE)


  G_TC <- localG(d_TC_poly$p, W_TC, zero.policy = TRUE)
  G_IC <- localG(d_IC_poly$p, W_IC, zero.policy = TRUE)



  localp_TC = attr(G_TC, "internal")[,5]; localp_TC[is.na(localp_TC)] = 1
  localp_IC = attr(G_IC, "internal")[,5]; localp_IC[is.na(localp_IC)] = 1



  ## both hotspot
  HS <- (attr(G_IC, "cluster") == "High" & localp_IC <= 0.05) & (attr(G_TC, "cluster") == "High" & localp_TC <= 0.05)
  ## both coldspot
  CS <- (attr(G_IC, "cluster") == "Low" & localp_IC <= 0.05) & (attr(G_TC, "cluster") == "Low" & localp_TC <= 0.05)
  ## tumor hotspot but immune coldspot
  CS_IC_HS_TC <- (attr(G_IC, "cluster") == "Low" & localp_IC <= 0.05) & (attr(G_TC, "cluster") == "High" & localp_TC <= 0.05)
  ## tumor coldspot but immune hotspot
  CS_TC_HS_IC <- (attr(G_IC, "cluster") == "High" & localp_IC <= 0.05) & (attr(G_TC, "cluster") == "Low" & localp_TC <= 0.05)

  ## immune hotspot
  HS_IC <- attr(G_IC, "cluster") == "High" & localp_IC <= 0.05
  ## immune coldspot
  CS_IC <- attr(G_IC, "cluster") == "Low" & localp_IC <= 0.05
  ## tumor hotspot
  HS_TC <- attr(G_TC, "cluster") == "High" & localp_TC <= 0.05
  ## tumor coldspot
  CS_TC <- attr(G_TC, "cluster") == "Low" & localp_TC <= 0.05


  GetisOrd_HS = sum(HS)/length(HS)
  GetisOrd_CS = sum(CS)/length(CS)
  GetisOrd_CS_IC_HS_TC = sum(CS_IC_HS_TC)/length(CS_IC_HS_TC)
  GetisOrd_CS_TC_HS_IC = sum(CS_TC_HS_IC)/length(CS_TC_HS_IC)

  GetisOrd_HS_IC = sum(HS_IC)/length(HS_IC)
  GetisOrd_CS_IC = sum(CS_IC)/length(CS_IC)
  GetisOrd_HS_TC = sum(HS_TC)/length(HS_TC)
  GetisOrd_CS_TC = sum(CS_TC)/length(CS_TC)

  GetisOrd_S_intra_immune = GetisOrd_HS/(GetisOrd_HS + GetisOrd_CS_TC_HS_IC + 0.000001) # add a small number to avoid 0 in the denominator
  GetisOrd_S_intra_cancer = GetisOrd_HS/(GetisOrd_HS + GetisOrd_CS_IC_HS_TC + 0.000001) # add a small number to avoid 0 in the denominator

  
  if (myplot){
    
    plot(d_TC_poly@data[, c("x", "y")], col=ifelse(HS_TC, "red", "black"), main = paste0(to_type, " hotspots"), axes=FALSE, xlab = "", ylab = "", frame.plot=F)
    plot(d_TC_poly@data[, c("x", "y")], col=ifelse(HS_IC, "red", "black"), main = paste0(from_type, " hotspots"), axes=FALSE, xlab = "", ylab = "", frame.plot=F)
    plot(d_TC_poly@data[, c("x", "y")], col=ifelse(HS, "red", "black"), main = paste0(to_type, " and ", from_type, " Co-localized hotspots"), axes=FALSE, xlab = "", ylab = "", frame.plot=F)

    plot(d_TC_poly@data[, c("x", "y")], col=ifelse(CS_TC, "blue", "black"), main = paste0(to_type, " coldspots"), axes=FALSE, xlab = "", ylab = "", frame.plot=F)
    plot(d_TC_poly@data[, c("x", "y")], col=ifelse(CS_IC, "blue", "black"), main = paste0(from_type, " coldspots"), axes=FALSE, xlab = "", ylab = "", frame.plot=F)
    plot(d_TC_poly@data[, c("x", "y")], col=ifelse(CS, "blue", "black"), main = paste0(to_type, " and ", from_type, " Co-localized coldspots"), axes=FALSE, xlab = "", ylab = "", frame.plot=F)

    plot(d_TC_poly@data[, c("x", "y")], col=ifelse(CS_IC_HS_TC, "red", "black"), main = paste0(to_type, " hotspots but ", from_type, " coldspots"), axes=FALSE, xlab = "", ylab = "", frame.plot=F)
    plot(d_TC_poly@data[, c("x", "y")], col=ifelse(CS_TC_HS_IC, "red", "black"), main = paste0(to_type, " coldspots but ", from_type, " hotspots"), axes=FALSE, xlab = "", ylab = "", frame.plot=F)
    
    
  }
  

  return(list("BC" = BC, "MH_index" = mh, "JaccardJ" = JJ, "SorensenL" = SL,
              "Moran_I_tumor" = moran_TC, "Moran_I_immune" = moran_IC, "moran_I_Bivariate" = moran_BV,
              "geary_TC" = geary_TC, "geary_IC" = geary_IC,
              "moran_HL_TC" = moran_HL_TC, "moran_HH_TC" = moran_HH_TC, "moran_LH_TC" = moran_LH_TC, "moran_LL_TC" = moran_LL_TC,
              "geary_HH_TC" = geary_HH_TC, "geary_LL_TC" = geary_LL_TC,
              "moran_HL_IC" = moran_HL_IC, "moran_HH_IC" = moran_HH_IC, "moran_LH_IC" = moran_LH_IC, "moran_LL_IC" = moran_LL_IC,
              "geary_HH_IC" = geary_HH_IC, "geary_LL_IC" = geary_LL_IC,
              "GetisOrd_HS" = GetisOrd_HS, "GetisOrd_CS" = GetisOrd_CS,
              "GetisOrd_CS_IC_HS_TC" = GetisOrd_CS_IC_HS_TC, "GetisOrd_CS_TC_HS_IC" = GetisOrd_CS_TC_HS_IC,
              "GetisOrd_HS_IC" = GetisOrd_HS_IC, "GetisOrd_CS_IC" = GetisOrd_CS_IC,
              "GetisOrd_HS_TC" = GetisOrd_HS_TC, "GetisOrd_CS_TC" = GetisOrd_CS_TC,
              "GetisOrd_S_intra_cancer" = GetisOrd_S_intra_cancer, "GetisOrd_S_intra_immune" = GetisOrd_S_intra_immune
  ))
}

