#' Main function to generate spatial features using point process data (bi-variate markers)
#'
#' This is a function that calculates the features for bi-variate markers, using spatial point process methods. 
#' The methods included are cross-type Ripley's K-function, differences in Ripley's K-function, cross-type G-function, Pair correlation function, Mark correlation function, Mark connection function, Marcon and Puech's M function.
#'
#' @param path The path for the directory that contains the data
#' @param file The file name in the directory
#' @param cell_class_var The column name in each file that indicates the cell class variable
#' @param x_var The column name in each file that indicates the cell location x variable
#' @param y_var The column name in each file that indicates the cell location y variable
#' @param from_type The cell type one wants to use as the "from" cell type
#' @param to_type The cell type one wants to use as the "to" cell type
#' @param scale The spatial range that user wants to investigate
#' @param myplot Whether to plot the results, if available, by default it is set as FALSE
#'
#' @return This function returns the features for bi-variate markers (denote from_type as type 1 and to_type as type 2), using spatial point process methods:
#' \item{g_cross_AUC}{AUC of cross-type G-function of type 1&2}
#' \item{g_cross_r}{the radius at which cross-type G-function reaches maximum, i.e. 1}
#' \item{k_cross_AUC}{AUC of cross-type Ripley's K-function of type 1&2}
#' \item{k_cross_vals_med}{median of cross-type K-function at specified range of radius}
#' \item{k_cross_vals_q1}{1st quantile of cross-type K-function at specified range of radius}
#' \item{k_cross_vals_q3}{3rd quantile of cross-type K-function at specified range of radius}
#' \item{k_cross_vals_max}{maximum of cross-type K-function at specified range of radius}
#' \item{K12_Diff_AUC}{difference between AUC of Ripley's K-function of type 1 and type 2}
#' \item{k_k1k2_vals_med}{median of difference between AUC of Ripley's K-function of type 1 and type 2 at specified range of radius}
#' \item{k_k1k2_vals_q1}{1st quantile of difference between AUC of Ripley's K-function of type 1 and type 2 at specified range of radius}
#' \item{k_k1k2_vals_q3}{3rd quantile of difference between AUC of Ripley's K-function of type 1 and type 2 at specified range of radius}
#' \item{k_k1k2_vals_max}{maximum of difference between AUC of Ripley's K-function of type 1 and type 2 at specified range of radius}
#' \item{K1K12_Diff_AUC}{difference between AUC of Ripley's K-function of type 1 and type 1&2}
#' \item{k_k1k12_vals_med}{median of difference between AUC of Ripley's K-function of type 1 and type 1&2 at specified range of radius}
#' \item{k_k1k12_vals_q1}{1st quantile of difference between AUC of Ripley's K-function of type 1 and type 1&2 at specified range of radius}
#' \item{k_k1k12_vals_q3}{3rd quantile of difference between AUC of Ripley's K-function of type 1 and type 1&2 at specified range of radius}
#' \item{k_k1k12_vals_max}{maximum of difference between AUC of Ripley's K-function of type 1 and type 1&2 at specified range of radius}
#' \item{K2K12_Diff_AUC}{difference between AUC of Ripley's K-function of type 2 and type 1&2}
#' \item{k_k2k12_vals_med}{median of difference between AUC of Ripley's K-function of type 2 and type 1&2 at specified range of radius}
#' \item{k_k2k12_vals_q1}{1st quantile of difference between AUC of Ripley's K-function of type 2 and type 1&2 at specified range of radius}
#' \item{k_k2k12_vals_q3}{3rd quantile of difference between AUC of Ripley's K-function of type 2 and type 1&2 at specified range of radius}
#' \item{k_k2k12_vals_max}{maximum of difference between AUC of Ripley's K-function of type 2 and type 1&2 at specified range of radius}
#' \item{pcf_AUC}{AUC of pair correlation function of type 1&2}
#' \item{pcf_vals_med}{median of pair correlation function at specified range of radius }
#' \item{pcf_vals_q1}{1st quantile of pair correlation function at specified range of radius }
#' \item{pcf_vals_q3}{3rd quantile of pair correlation function at specified range of radius }
#' \item{pcf_vals_max}{maximum of pair correlation function at specified range of radius }
#' \item{pcf_vals_min}{minimum of pair correlation function at specified range of radius }
#' \item{pcf_r}{the radius at which pair correlation function reaches maximum}
#' \item{mrcf_AUC}{AUC of mark correlation function of type 1&2}
#' \item{mrcf_vals_med}{median of mark correlation function at specified range of radius }
#' \item{mrcf_vals_q1}{1st quantile of mark correlation function at specified range of radius }
#' \item{mrcf_vals_q3}{3rd quantile of mark correlation function at specified range of radius }
#' \item{mrcf_vals_max}{maximum of mark correlation function at specified range of radius }
#' \item{mrcf_vals_min}{minimum of mark correlation function at specified range of radius }
#' \item{mrcf_r}{the radius at which mark correlation function reaches maximum}
#' \item{mccf_AUC}{AUC of mark connection function of type 1&2}
#' \item{mccf_vals_med}{median of mark connection function at specified range of radius }
#' \item{mccf_vals_q1}{1st quantile of mark connection function at specified range of radius }
#' \item{mccf_vals_q3}{3rd quantile of mark connection function at specified range of radius }
#' \item{mccf_vals_max}{maximum of mark connection function at specified range of radius }
#' \item{mccf_vals_min}{minimum of mark connection function at specified range of radius }
#' \item{mccf_r}{the radius at which mark connection function reaches maximum}
#' \item{M_AUC}{AUC of Marcon and Puech'S M function of type 1&2}
#' \item{M_vals_med}{median of Marcon and Puech'S M function at specified range of radius }
#' \item{M_vals_q1}{1st quantile of Marcon and Puech'S M function at specified range of radius }
#' \item{M_vals_q3}{3rd quantile of Marcon and Puech'S M function at specified range of radius }
#' \item{M_vals_max}{maximum of Marcon and Puech'S M function at specified range of radius }
#' \item{M_vals_min}{minimum of Marcon and Puech'S M function at specified range of radius }
#' \item{M_r}{the radius at which Marcon and Puech'S M function reaches maximum}
#'
#' @import magrittr
#' @import dplyr
#' @import spatstat.geom
#' @import spatstat.explore
#' @import spatstat.utils
#' @import ggplot2
#' @import rpart
#' @importFrom dbmss wmppp Mhat
#' @importFrom ecespa K1K2
#' @importFrom pracma trapz
#' @importFrom sm sm.density
#'
#' @author Xiao Li, \email{xiao.li.xl2@roche.com}
#'
#' @export







Point_pattern_data_bi <- function(path = "/Users/lix233/Haystack/5862_cell_centers/",
                                       file = "f8008cf4-3fe7-4200-a2af-21c27708ae28.csv",
                                       cell_class_var = "cell_class",
                                       x_var = "X0",
                                       y_var = "X1",
                                       from_type = "Lymphocyte",
                                       to_type = "Tumor",
                                       scale,
                                       myplot = FALSE){


  spp_df = read.csv(paste0(path, file))

  spp_df = spp_df %>% select("cell_class" = all_of(cell_class_var),
                             "x_coord" = all_of(x_var),
                             "y_coord" = all_of(y_var))

  L = max(max(spp_df$x_coord), max(spp_df$y_coord))

  spp_df = spp_df %>% mutate("x" = normalized_coords(spp_df$x_coord, L),
                             "y" = normalized_coords(spp_df$y_coord, L),
                             "cell_class" = as.factor(spp_df$cell_class))


  radius = scale/L


  ##########################################################################################################################################
  ##########################################################################################################################################
  ####################################################        Spatial plots & features        ##############################################
  ##########################################################################################################################################
  ##########################################################################################################################################

  # generate the enclosing polygon
  xy <- as.matrix(spp_df[, c('x','y')])

  bnd = owin(xrange = c(0, max(xy[,1])), yrange = c(0, max(xy[,2])))


  ln = with(spp_df,
            ppp(x = x, y = y, marks = cell_class, window = bnd)
  )


  spp_immune = spp_df %>% filter(cell_class == from_type)
  spp_tumor = spp_df %>% filter(cell_class == to_type)


  ## Cross type G-function
  Gln = Gcross(ln, i = from_type, j = to_type, r=seq(0, radius, length.out = 10), correction = "km")

  if (myplot) {
    plot(Gln, lty = 1, main = paste0("Corss type G-function for ", from_type, " and ", to_type))
  }   
  tmp = as.data.frame(cbind(Gln$km, Gln$r, Gln$theo)); colnames(tmp) = c('km', 'r', 'theo')
  tmp = tmp[is.finite(tmp$km), ]
  g_cross_AUC = pracma::trapz(tmp$r, tmp$km-tmp$theo)
  g_cross_r = tmp$r[which.max(tmp$km)]
  
  cat("Cross type G-function is calculated. \n")


  ## Corss type K-function
  Kln = Kcross(ln, i = from_type, j = to_type, r=seq(0, radius, length.out = 10), correction = "translate")

  if (myplot) {
    plot(Kln, lty = 1, main = paste0("Corss type K-function for ", from_type, " and ", to_type))
  } 
  tmp = as.data.frame(cbind(Kln$trans, Kln$r, Kln$theo)); colnames(tmp) = c('trans', 'r', 'theo')
  tmp = tmp[is.finite(tmp$trans), ]
  k_cross_AUC = pracma::trapz(tmp$r, tmp$trans-tmp$theo)
  k_cross_vals_med = quantile(tmp$trans, 0.5)
  k_cross_vals_q1 = quantile(tmp$trans, 0.25)
  k_cross_vals_q3 = quantile(tmp$trans, 0.75)
  k_cross_vals_max = max(tmp$trans)
  
  
  cat("Cross type K-function is calculated. \n")



  ## differences between univariate and bivariate K-functions
  K_diff = ecespa::K1K2(ln, i = from_type, j = to_type, nsim = 2, nrank = 1, r=seq(0, radius, length.out = 10),
                correction = "translate")
  
  if (myplot) {
    plot(K_diff$k1k2, lty=c(2, 1, 2), col=c(2, 1, 2), main = paste0("Difference between K-function for ", from_type, " and ", to_type))
  } 

  # tmp = as.data.frame(cbind(K_diff$k1k2$`K1-K2`, K_diff$k1k2$r)); colnames(tmp) = c('trans', 'r')
  tmp = as.data.frame(cbind(K_diff$k1k2$D, K_diff$k1k2$r)); colnames(tmp) = c('trans', 'r')
  tmp = tmp[is.finite(tmp$trans), ]
  k_k1k2_AUC = pracma::trapz(tmp$r, tmp$trans)
  k_k1k2_vals_med = quantile(tmp$trans, 0.5)
  k_k1k2_vals_q1 = quantile(tmp$trans, 0.25)
  k_k1k2_vals_q3 = quantile(tmp$trans, 0.75)
  k_k1k2_vals_max = max(tmp$trans)

  if (myplot) {
    plot(K_diff$k1k12, lty=c(2, 1, 2), col=c(2, 1, 2), main = paste0("Difference between K-function for ", from_type, " and ", to_type, " vs.", from_type))
  } 
  tmp = as.data.frame(cbind(K_diff$k1k12$D, K_diff$k1k12$r)); colnames(tmp) = c('trans', 'r')
  tmp = tmp[is.finite(tmp$trans), ]
  k_k1k12_AUC = pracma::trapz(tmp$r, tmp$trans)
  k_k1k12_vals_med = quantile(tmp$trans, 0.5)
  k_k1k12_vals_q1 = quantile(tmp$trans, 0.25)
  k_k1k12_vals_q3 = quantile(tmp$trans, 0.75)
  k_k1k12_vals_max = max(tmp$trans)

  if (myplot) {
    plot(K_diff$k2k12, lty=c(2, 1, 2), col=c(2, 1, 2), main = paste0("Difference between K-function for ", to_type, " and ", to_type, " vs.", from_type))
  } 
  tmp = as.data.frame(cbind(K_diff$k2k12$D, K_diff$k1k12$r)); colnames(tmp) = c('trans', 'r')
  tmp = tmp[is.finite(tmp$trans), ]
  k_k2k12_AUC = pracma::trapz(tmp$r, tmp$trans)
  k_k2k12_vals_med = quantile(tmp$trans, 0.5)
  k_k2k12_vals_q1 = quantile(tmp$trans, 0.25)
  k_k2k12_vals_q3 = quantile(tmp$trans, 0.75)
  k_k2k12_vals_max = max(tmp$trans)
  
  
  cat("Differences in Ripley's K-function is calculated. \n")




  pcfmulti <-
    function (X, I, J, ..., r = NULL, kernel = "epanechnikov", bw = NULL,
              stoyan = 0.15, correction = c("translate", "Ripley"), divisor = c("r",
                                                                                "d"), Iname = "points satisfying condition I", Jname = "points satisfying condition J")
    {
      verifyclass(X, "ppp")
      divisor <- match.arg(divisor)
      win <- X$window
      areaW <- spatstat.geom::area(win)
      npts <- npoints(X)
      correction.given <- !missing(correction) && !is.null(correction)
      if (is.null(correction))
        correction <- c("translate", "Ripley")
      correction <- pickoption("correction", correction, c(isotropic = "isotropic",
                                                           Ripley = "isotropic", trans = "translate", translate = "translate",
                                                           translation = "translate", best = "best"), multi = TRUE)
      correction <- implemented.for.K(correction, win$type, correction.given)
      I <- ppsubset(X, I)
      J <- ppsubset(X, J)
      if (is.null(I) || is.null(J))
        stop("I and J must be valid subset indices")
      nI <- sum(I)
      nJ <- sum(J)
      if (nI == 0)
        stop(paste("There are no", Iname))
      if (nJ == 0)
        stop(paste("There are no", Jname))
      XI <- X[I]
      XJ <- X[J]
      lambdaJ <- nJ/areaW
      nIJ <- sum(I & J)
      lambdaIJarea <- (as.numeric(nI) * as.numeric(nJ) - as.numeric(nIJ))/areaW
      if (is.null(bw) && kernel == "epanechnikov") {
        h <- stoyan/sqrt(lambdaJ)
        hmax <- h
        bw <- h/sqrt(5)
      }
      else if (is.numeric(bw)) {
        hmax <- 3 * bw
      }
      else {
        hmax <- 2 * stoyan/sqrt(lambdaJ)
      }
      rmaxdefault <- rmax.rule("K", win, lambdaJ)
      breaks <- handle.r.b.args(r, NULL, win, rmaxdefault = rmaxdefault)
      if (!(breaks$even))
        stop("r values must be evenly spaced")
      r <- breaks$r
      rmax <- breaks$max
      alim <- c(0, min(rmax, rmaxdefault))
      df <- data.frame(r = r, theo = rep.int(1, length(r)))
      fname <- c("g", "list(I,J)")
      yexp <- quote(g[list(I, J)](r))
      out <- fv(df, "r", quote(g[I, J](r)), "theo", , alim, c("r",
                                                              makefvlabel(NULL, NULL, fname, "Pois")), c("distance argument r",
                                                                                                         "theoretical Poisson %s"), fname = fname, yexp = yexp)
      denargs <- spatstat.utils::resolve.defaults(list(kernel = kernel, bw = bw),
                                                  list(...), list(n = length(r), from = 0, to = rmax))
      what <- if (any(correction == "translate"))
        "all"
      else "ijd"
      close <- crosspairs(XI, XJ, rmax + hmax, what = what)
      orig <- seq_len(npts)
      imap <- orig[I]
      jmap <- orig[J]
      iX <- imap[close$i]
      jX <- jmap[close$j]
      if (nIJ > 0) {
        ok <- (iX != jX)
        if (!all(ok))
          close <- as.list(as.data.frame(close)[ok, , drop = FALSE])
      }
      dclose <- close$d
      icloseI <- close$i
      if (any(correction == "translate")) {
        edgewt <- edge.Trans(dx = close$dx, dy = close$dy, W = win,
                             paired = TRUE)
        gT <- sewpcf(dclose, edgewt, denargs, lambdaIJarea, divisor)$g
        out <- bind.fv(out, data.frame(trans = gT), makefvlabel(NULL,
                                                                "hat", fname, "Trans"), "translation-corrected estimate of %s",
                       "trans")
      }
      if (any(correction == "isotropic")) {
        edgewt <- edge.Ripley(XI[icloseI], matrix(dclose, ncol = 1))
        gR <- sewpcf(dclose, edgewt, denargs, lambdaIJarea, divisor)$g
        out <- bind.fv(out, data.frame(iso = gR), makefvlabel(NULL,
                                                              "hat", fname, "Ripley"), "isotropic-corrected estimate of %s",
                       "iso")
      }
      if (is.null(out)) {
        warning("Nothing computed - no edge corrections chosen")
        return(NULL)
      }
      corrxns <- rev(setdiff(names(out), "r"))
      formula(out) <- . ~ r
      fvnames(out, ".") <- corrxns
      unitname(out) <- unitname(X)
      return(out)
    }


  pcfcross <-
    function (X, i, j, ..., r = NULL, kernel = "epanechnikov", bw = NULL,
              stoyan = 0.15, correction = c("isotropic", "Ripley", "translate"),
              divisor = c("r", "d"))
    {
      verifyclass(X, "ppp")
      stopifnot(is.multitype(X))
      if (missing(correction))
        correction <- NULL
      divisor <- match.arg(divisor)
      marx <- marks(X)
      if (missing(i))
        i <- levels(marx)[1]
      if (missing(j))
        j <- levels(marx)[2]
      I <- (marx == i)
      J <- (marx == j)
      Iname <- paste("points with mark i =", i)
      Jname <- paste("points with mark j =", j)
      result <- pcfmulti(X, I, J, ..., r = r, kernel = kernel,
                         bw = bw, stoyan = stoyan, correction = correction, divisor = divisor,
                         Iname = Iname, Jname = Jname)
      iname <- spatstat.utils::make.parseable(paste(i))
      jname <- spatstat.utils::make.parseable(paste(j))
      result <- rebadge.fv(result, substitute(g[i, j](r), list(i = iname,
                                                               j = jname)), c("g", paste0("list", spatstat.utils::paren(paste(iname,
                                                                                                                              jname, sep = ",")))), new.yexp = substitute(g[list(i,
                                                                                                                                                                                 j)](r), list(i = iname, j = jname)))
      return(result)
    }




  sewsmod <- function(d, ff, wt, Ef, rvals, method="smrep", ..., nwtsteps=500) {
    ## Smooth Estimate of Weighted Second Moment Density
    ## (engine for computing mark correlations, etc)
    ## ------
    ## Vectors containing one entry for each (close) pair of points
    ## d = interpoint distance
    ## ff = f(M1, M2) where M1, M2 are marks at the two points
    ## wt = edge correction weight
    ## -----
    ## Ef = E[f(M, M')] where M, M' are independent random marks
    ##
    d <- as.vector(d)
    ff <- as.vector(ff)
    wt <- as.vector(wt)
    switch(method,
           density={
             fw <- ff * wt
             sum.fw <- sum(fw)
             sum.wt <- sum(wt)
             ## smooth estimate of kappa_f
             est <- density(d, weights=fw/sum.fw,
                            from=min(rvals), to=max(rvals), n=length(rvals),
                            ...)$y
             numerator <- est * sum.fw
             ## smooth estimate of kappa_1
             est0 <- density(d, weights=wt/sum.wt,
                             from=min(rvals), to=max(rvals), n=length(rvals),
                             ...)$y
             denominator <- est0 * Ef * sum.wt
             result <- numerator/denominator
           },
           sm={
             ## This is slow!
             suppressWarnings(smok <- requireNamespace("sm"))
             if(!smok)
               stop(paste("Option method=sm requires package sm,",
                          "which is not available"))

             ## smooth estimate of kappa_f
             fw <- ff * wt
             est <- sm::sm.density(d, weights=fw,
                                   eval.points=rvals,
                                   display="none", nbins=0, ...)$estimate
             numerator <- est * sum(fw)/sum(est)
             ## smooth estimate of kappa_1
             est0 <- sm::sm.density(d, weights=wt,
                                    eval.points=rvals,
                                    display="none", nbins=0, ...)$estimate
             denominator <- est0 * (sum(wt)/ sum(est0)) * Ef
             result <- numerator/denominator
           },
           smrep={
             suppressWarnings(smok <- requireNamespace("sm"))
             if(!smok)
               stop(paste("Option method=smrep requires package sm,",
                          "which is not available"))

             hstuff <- resolve.defaults(list(...), list(hmult=1, h.weights=NA))
             if(hstuff$hmult == 1 && all(is.na(hstuff$h.weights)))
               warning("default smoothing parameter may be inappropriate")

             ## use replication to effect the weights (it's faster)
             nw <- round(nwtsteps * wt/max(wt))
             drep.w <- rep.int(d, nw)
             fw <- ff * wt
             nfw <- round(nwtsteps * fw/max(fw))
             drep.fw <- rep.int(d, nfw)

             ## smooth estimate of kappa_f
             est <- sm::sm.density(drep.fw,
                                   eval.points=rvals,
                                   display="none", ...)$estimate
             numerator <- est * sum(fw)/sum(est)
             ## smooth estimate of kappa_1
             est0 <- sm::sm.density(drep.w,
                                    eval.points=rvals,
                                    display="none", ...)$estimate
             denominator <- est0 * (sum(wt)/ sum(est0)) * Ef
             result <- numerator/denominator
           },
           loess = {
             ## set up data frame
             df <- data.frame(d=d, ff=ff, wt=wt)
             ## fit curve to numerator using loess
             fitobj <- loess(ff ~ d, data=df, weights=wt, ...)
             ## evaluate fitted curve at desired r values
             Eff <- predict(fitobj, newdata=data.frame(d=rvals))
             ## normalise:
             ## denominator is the sample mean of all ff[i,j],
             ## an estimate of E(ff(M1,M2)) for M1,M2 independent marks
             result <- Eff/Ef
           },
    )
    return(result)
  }


  markcorr <-
    function (X, f = function(m1, m2) {
      m1 * m2
    }, r = NULL, correction = c("isotropic", "Ripley", "translate"),
    method = "density", ..., weights = NULL, f1 = NULL, normalise = TRUE,
    fargs = NULL, internal = NULL)
    {
      stopifnot(is.ppp(X) && is.marked(X))
      nX <- npoints(X)
      if (missing(f))
        f <- NULL
      if (missing(correction))
        correction <- NULL
      marx <- marks(X, dfok = TRUE)
      if (is.data.frame(marx)) {
        nc <- ncol(marx)
        result <- list()
        for (j in 1:nc) {
          Xj <- X %mark% marx[, j]
          result[[j]] <- markcorr(Xj, f = f, r = r, correction = correction,
                                  method = method, ..., weights = weights, f1 = f1,
                                  normalise = normalise, fargs = fargs)
        }
        result <- as.anylist(result)
        names(result) <- colnames(marx)
        return(result)
      }
      if (unweighted <- is.null(weights)) {
        weights <- rep(1, nX)
      }
      else {
        weights <- pointweights(X, weights = weights, parent = parent.frame())
        stopifnot(all(weights > 0))
      }
      h <- check.testfun(f, f1, X)
      f <- h$f
      f1 <- h$f1
      ftype <- h$ftype
      npts <- npoints(X)
      W <- X$window
      rmaxdefault <- rmax.rule("K", W, npts/area(W))
      breaks <- handle.r.b.args(r, NULL, W, rmaxdefault = rmaxdefault)
      r <- breaks$r
      rmax <- breaks$max
      if (length(method) > 1)
        stop("Select only one method, please")
      if (method == "density" && !breaks$even)
        stop(paste("Evenly spaced r values are required if method=",
                   sQuote("density"), sep = ""))
      correction.given <- !missing(correction) && !is.null(correction)
      if (is.null(correction))
        correction <- c("isotropic", "Ripley", "translate")
      correction <- pickoption("correction", correction, c(none = "none",
                                                           border = "border", bord.modif = "bord.modif", isotropic = "isotropic",
                                                           Ripley = "isotropic", translate = "translate", translation = "translate",
                                                           best = "best"), multi = TRUE)
      correction <- implemented.for.K(correction, W$type, correction.given)
      Ef <- internal$Ef
      if (is.null(Ef)) {
        Ef <- switch(ftype, mul = {
          mean(marx * weights)^2
        }, equ = {
          if (unweighted) {
            mtable <- table(marx)
          } else {
            mtable <- tapply(weights, marx, sum)
            mtable[is.na(mtable)] <- 0
          }
          sum(mtable^2)/nX^2
        }, product = {
          f1m <- do.call(f1, append(list(marx), fargs))
          mean(f1m * weights)^2
        }, general = {
          mcross <- if (is.null(fargs)) {
            outer(marx, marx, f)
          } else {
            do.call(outer, append(list(marx, marx, f), fargs))
          }
          if (unweighted) {
            mean(mcross)
          } else {
            wcross <- outer(weights, weights, "*")
            mean(mcross * wcross)
          }
        }, stop("Internal error: invalid ftype"))
      }
      if (normalise) {
        theory <- 1
        Efdenom <- Ef
      }
      else {
        theory <- Ef
        Efdenom <- 1
      }
      if (normalise) {
        if (Efdenom == 0)
          stop("Cannot normalise the mark correlation; the denominator is zero")
        else if (Efdenom < 0)
          warning(paste("Problem when normalising the mark correlation:",
                        "the denominator is negative"))
      }
      result <- data.frame(r = r, theo = rep.int(theory, length(r)))
      desc <- c("distance argument r", "theoretical value (independent marks) for %s")
      alim <- c(0, min(rmax, rmaxdefault))
      if (ftype %in% c("mul", "equ")) {
        if (normalise) {
          ylab <- quote(k[mm](r))
          fnam <- c("k", "mm")
        }
        else {
          ylab <- quote(c[mm](r))
          fnam <- c("c", "mm")
        }
      }
      else {
        if (normalise) {
          ylab <- quote(k[f](r))
          fnam <- c("k", "f")
        }
        else {
          ylab <- quote(c[f](r))
          fnam <- c("c", "f")
        }
      }
      result <- fv(result, "r", ylab, "theo", , alim, c("r", "{%s[%s]^{iid}}(r)"),
                   desc, fname = fnam)
      close <- closepairs(X, rmax)
      dIJ <- close$d
      I <- close$i
      J <- close$j
      XI <- ppp(close$xi, close$yi, window = W, check = FALSE)
      mI <- marx[I]
      mJ <- marx[J]
      ff <- switch(ftype, mul = mI * mJ, equ = (mI == mJ), product = {
        if (is.null(fargs)) {
          fI <- f1(mI)
          fJ <- f1(mJ)
        } else {
          fI <- do.call(f1, append(list(mI), fargs))
          fJ <- do.call(f1, append(list(mJ), fargs))
        }
        fI * fJ
      }, general = {
        if (is.null(fargs)) f(marx[I], marx[J]) else do.call(f,
                                                             append(list(marx[I], marx[J]), fargs))
      })
      if (is.logical(ff))
        ff <- as.numeric(ff)
      else if (!is.numeric(ff))
        stop("function f did not return numeric values")
      if (anyNA(ff))
        switch(ftype, mul = , equ = stop("some marks were NA"),
               product = , general = stop("function f returned some NA values"))
      if (any(ff < 0))
        switch(ftype, mul = , equ = stop("negative marks are not permitted"),
               product = , general = stop("negative values of function f are not permitted"))
      if (!unweighted)
        ff <- ff * weights[I] * weights[J]
      if (any(correction == "none")) {
        edgewt <- rep.int(1, length(dIJ))
        Mnone <- sewsmod(dIJ, ff, edgewt, Efdenom, r, method,
                         ...)
        result <- bind.fv(result, data.frame(un = Mnone), "{hat(%s)[%s]^{un}}(r)",
                          "uncorrected estimate of %s", "un")
      }
      if (any(correction == "translate")) {
        XJ <- ppp(close$xj, close$yj, window = W, check = FALSE)
        edgewt <- edge.Trans(XI, XJ, paired = TRUE)
        Mtrans <- sewsmod(dIJ, ff, edgewt, Efdenom, r, method,
                          ...)
        result <- bind.fv(result, data.frame(trans = Mtrans),
                          "{hat(%s)[%s]^{trans}}(r)", "translation-corrected estimate of %s",
                          "trans")
      }
      if (any(correction == "isotropic")) {
        edgewt <- edge.Ripley(XI, matrix(dIJ, ncol = 1))
        Miso <- sewsmod(dIJ, ff, edgewt, Efdenom, r, method,
                        ...)
        result <- bind.fv(result, data.frame(iso = Miso), "{hat(%s)[%s]^{iso}}(r)",
                          "Ripley isotropic correction estimate of %s", "iso")
      }
      nama2 <- names(result)
      corrxns <- rev(nama2[nama2 != "r"])
      formula(result) <- (. ~ r)
      fvnames(result, ".") <- corrxns
      unitname(result) <- unitname(X)
      return(result)
    }



  markconnect <- local({

    indicateij <- function(m1, m2, i, j) { (m1 == i) & (m2 == j) }

    markconnect <- function(X, i, j, r=NULL,
                            correction=c("isotropic", "Ripley", "translate"),
                            method="density", ...,
                            normalise=FALSE) {
      stopifnot(is.ppp(X) && is.multitype(X))
      if(missing(correction))
        correction <- NULL
      marx <- marks(X)
      lev  <- levels(marx)
      if(missing(i)) i <- lev[1]
      if(missing(j)) j <- lev[2]
      ## compute reference value Ef
      weights <- pointweights(X, ..., parent=parent.frame())
      Ef <- if(is.null(weights)) mean(marx == i) * mean(marx == j) else
        mean(weights * (marx == i)) * mean(weights * (marx == j))
      ## compute estimates
      p <- markcorr(X, f=indicateij, r=r,
                    correction=correction, method=method,
                    ...,
                    fargs=list(i=i, j=j),
                    normalise=normalise,
                    internal=list(Ef=Ef))
      ## alter theoretical value and fix labels
      if(!normalise) {
        pipj <- mean(marx==i) * mean(marx==j)
        p$theo <- pipj
      } else {
        p$theo <- 1
      }
      p <- rebadge.fv(p,
                      new.ylab=substitute(p[i,j](r), list(i=paste(i),j=paste(j))),
                      new.fname=c("p", paste0("list(", i, ",", j, ")")),
                      new.yexp=substitute(p[list(i,j)](r),
                                          list(i=paste(i),j=paste(j))))
      return(p)
    }
    markconnect
  })









  ## cross type PCF
  pcfln = pcfcross(ln, r = seq(0, radius, length.out = 10), i = from_type, j = to_type, correction = "translate")

  if (myplot) {
    plot(pcfln, lty=1, main = paste0("Pair correlation function for ", from_type, " and ", to_type))
  }

  tmp = as.data.frame(cbind(pcfln$trans, pcfln$r, pcfln$theo)); colnames(tmp) = c('trans', 'r', 'theo')
  tmp = tmp[is.finite(tmp$trans), ]
  pcf_AUC = pracma::trapz(tmp$r, tmp$trans - tmp$theo)
  pcf_vals_med = quantile(tmp$trans, 0.5)
  pcf_vals_q1 = quantile(tmp$trans, 0.25)
  pcf_vals_q3 = quantile(tmp$trans, 0.75)
  pcf_vals_max = max(tmp$trans)
  pcf_vals_min = min(tmp$trans)
  pcf_r = tmp$r[which.max(tmp$trans)]
  
  cat("Pair correlation function is calculated. \n")

  ## mark correlation function
  mrcfln = markcorr(ln, r=seq(0, radius, length.out = 10), correction = "translate")

  if (myplot) {
    plot(mrcfln, lty=1, main = paste0("Mark correlation function for ", from_type, " and ", to_type))
  }
  tmp = as.data.frame(cbind(mrcfln$trans, mrcfln$r, mrcfln$theo)); colnames(tmp) = c('trans', 'r', 'theo')
  tmp = tmp[is.finite(tmp$trans), ]
  mrcf_AUC = pracma::trapz(tmp$r, tmp$trans - tmp$theo)
  mrcf_vals_med = quantile(tmp$trans, 0.5)
  mrcf_vals_q1 = quantile(tmp$trans, 0.25)
  mrcf_vals_q3 = quantile(tmp$trans, 0.75)
  mrcf_vals_max = max(tmp$trans)
  mrcf_vals_min = min(tmp$trans)
  mrcf_r = tmp$r[which.max(tmp$trans)]
  
  cat("Mark correlation function is calculated.")

  ## mark connection function
  mccfln = markconnect(ln, r=seq(0, radius, length.out = 10), i = from_type, j = to_type, correction = "translate")

  if (myplot) {
    plot(mccfln, lty=1, main = paste0("Mark connection function for ", from_type, " and ", to_type))
  }
  tmp = as.data.frame(cbind(mccfln$trans, mccfln$r, mccfln$theo)); colnames(tmp) = c('trans', 'r', 'theo')
  tmp = tmp[is.finite(tmp$trans), ]
  mccf_AUC = pracma::trapz(tmp$r, tmp$trans - tmp$theo)
  mccf_vals_med = quantile(tmp$trans, 0.5)
  mccf_vals_q1 = quantile(tmp$trans, 0.25)
  mccf_vals_q3 = quantile(tmp$trans, 0.75)
  mccf_vals_max = max(tmp$trans)
  mccf_vals_min = min(tmp$trans)
  mccf_r = tmp$r[which.max(tmp$trans)]
  
  cat("Mark connection function is calculated. \n")

  ## Marcon and Puech'S M function
  lnn = dbmss::wmppp(data.frame("X" = spp_df$x, "Y" = spp_df$y, "PointType" = spp_df$cell_class), window = bnd)
  Mfun = dbmss::Mhat(lnn, ReferenceType = from_type, NeighborType = to_type, CaseControl = TRUE, r=seq(0, radius, length.out = 10))

  if (myplot) {
    plot(Mfun, lty=1, main = paste0("Marcon and Puech'S M function for ", from_type, " and ", to_type))
  }
  tmp = as.data.frame(cbind(Mfun$M, Mfun$r, Mfun$theo)); colnames(tmp) = c('M', 'r', 'theo')
  tmp = tmp[is.finite(tmp$M), ]
  M_AUC = pracma::trapz(tmp$r, tmp$M - tmp$theo)
  M_vals_med = quantile(tmp$M, 0.5)
  M_vals_q1 = quantile(tmp$M, 0.25)
  M_vals_q3 = quantile(tmp$M, 0.75)
  M_vals_max = max(tmp$M)
  M_vals_min = min(tmp$M)
  M_r = tmp$r[which.max(tmp$M)]

  cat("Marcon and Puech's M function is calculated. \n")
  

  return(list("g_cross_AUC" = g_cross_AUC,
              "g_cross_r" = g_cross_r,
              
              "k_cross_AUC" = k_cross_AUC,
              "k_cross_vals_med" = k_cross_vals_med,
              "k_cross_vals_q1" = k_cross_vals_q1,
              "k_cross_vals_q3" = k_cross_vals_q3,
              "k_cross_vals_max" = k_cross_vals_max,
              
              "K12_Diff_AUC" = k_k1k2_AUC,
              "k_k1k2_vals_med" = k_k1k2_vals_med,
              "k_k1k2_vals_q1" = k_k1k2_vals_q1,
              "k_k1k2_vals_q3" = k_k1k2_vals_q3,
              "k_k1k2_vals_max" = k_k1k2_vals_max,
              
              "K1K12_Diff_AUC" = k_k1k12_AUC,
              "k_k1k12_vals_med" = k_k1k12_vals_med,
              "k_k1k12_vals_q1" = k_k1k12_vals_q1,
              "k_k1k12_vals_q3" = k_k1k12_vals_q3,
              "k_k1k12_vals_max" = k_k1k12_vals_max,
                            
              "K2K12_Diff_AUC" = k_k2k12_AUC,
              "k_k2k12_vals_med" = k_k2k12_vals_med,
              "k_k2k12_vals_q1" = k_k2k12_vals_q1,
              "k_k2k12_vals_q3" = k_k2k12_vals_q3,
              "k_k2k12_vals_max" = k_k2k12_vals_max,
              
              "pcf_AUC" = pcf_AUC,
              "pcf_vals_med" = pcf_vals_med,
              "pcf_vals_q1" = pcf_vals_q1,
              "pcf_vals_q3" = pcf_vals_q3,
              "pcf_vals_max" = pcf_vals_max,
              "pcf_vals_min" = pcf_vals_min,
              "pcf_r" = pcf_r,
              
              "mrcf_AUC" = mrcf_AUC,
              "mrcf_vals_med" = mrcf_vals_med,
              "mrcf_vals_q1" = mrcf_vals_q1,
              "mrcf_vals_q3" = mrcf_vals_q3,
              "mrcf_vals_max" = mrcf_vals_max,
              "mrcf_vals_min" = mrcf_vals_min,
              "mrcf_r" = mrcf_r,
              
              "mccf_AUC" = mccf_AUC,
              "mccf_vals_med" = mccf_vals_med,
              "mccf_vals_q1" = mccf_vals_q1,
              "mccf_vals_q3" = mccf_vals_q3,
              "mccf_vals_max" = mccf_vals_max,
              "mccf_vals_min" = mccf_vals_min,
              "mccf_r" = mccf_r,
              
              "M_AUC" = M_AUC,
              "M_vals_med" = M_vals_med,
              "M_vals_q1" = M_vals_q1,
              "M_vals_q3" = M_vals_q3,
              "M_vals_max" = M_vals_max,
              "M_vals_min" = M_vals_min,
              "M_r" = M_r))

}
