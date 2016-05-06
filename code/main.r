# Saharan Dust Transport Research
# 2016-03-19 Ismail SEZEN
# sezenismail@gmail.com
# levels: 1000  925  850  700  600  500  400  300  250
#          200  150  100   70   50   30   20   10

library(pryr)
library(tools)
library(rgdal)
library(rpart)
library(rpart.plot)
library(plyr)
library(ggplot2)
library(xts)
# library(lubridate)

library(plotrix)
library(ncdf4)
library(ncdf4.helpers)
library(TTR)
library(timeSeries)

library(akima)
library(rasterVis)

source("code/filehelper.r")

# INITIALS
dir_data <- "data"
dir_out  <- "out"

normalize <- function(x, perct=1, in_range=c(0, 1), remove_mean=F) {
  xs <- as.numeric(sort(x))
  i <- round(perct * length(xs))
  mx <- as.numeric(xs[i])
  x[x > mx] <- mx
  r <- range(x, na.rm = T)
  X <- (x - r[1]) / (r[2] - r[1])
  Xs <- X * (max(in_range) - min(in_range)) + min(in_range)
  if (remove_mean) Xs <- Xs - mean(Xs, na.rm = T)
  return(Xs)
}

create_grids <- function(x, y, dx=NULL, dy = dx) {
  if (is.null(dx)) dx <- x[2] - x[1]
  if (is.null(dy)) dy <- y[2] - y[1]
  xx <- seq(from = max(x), to = min(x), by = -abs(dx))
  yy <- seq(from = min(y), to = max(y), by = abs(dy))
  lx <- length(xx); ly <- length(yy)
  xx <- rep(xx, ly)
  yy <- rep(yy, each = lx)
  res <- cbind(xx, yy)[len(xx):1, ]
  cn <- c(deparse(substitute(x)), deparse(substitute(y)))
  colnames(res) <- cn
  return(res)
}

plot_streamlines <- function(u, v, ...) {
  l <- interpolate( u, v)
  u <- l$u; v <- l$v
  lon <- as.numeric(rownames(u))
  lat <- as.numeric(colnames(u))
  xmn = min(lon); xmx = max(lon)
  ymn = min(lat); ymx = max(lat)
  u1 <- raster(x = t(u), xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx)
  v1 <- raster(x = t(v), xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx)
  w <- stack(u1, v1)
  p <- streamplot(w, isField = 'dXY', region = F, droplet = list(pc = .3),
             xlab = "Longitude", ylab = "Latitude", lwd = 1.2,
             par.settings = streamTheme(
               panel.background = list(col = 'cornflowerblue'),
               symbol = brewer.pal(n = 9, name = 'Reds')),
             streamlet = list(L = 30), ...) +
    latticeExtra::layer(sp.polygons(cntry, fill = 'gainsboro',
                                    col = "black", alpha = 0.8), under = T)
  return(p)
}

stackplot <- function(data, ylim=NA, main=NA, colors=NA, xlab=NA, ylab=NA) {
  # stacked line plot
  if (is.na(ylim)) {
    ylim <- c(-50, max(rowSums(data, na.rm = T)))
  }
  if (is.na(colors)) {
    colors <- c("green", "red", "lightgray", "blue",
               "orange", "purple", "yellow")
  }
  xval <- as.numeric(row.names(data))
  summary <- rep(0, nrow(data))
  recent <- summary

  # Create empty plot
  plot(c(-100), c(-100), xlim = c(min(xval, na.rm = T), max(xval, na.rm = T)),
       ylim = ylim, main = main, xlab = xlab, ylab = ylab)

  # One polygon per column
  cols <- names(data)
  for (c in 1:length(cols)) {
    current <- data[[cols[[c]]]]
    summary <- summary + current
    polygon(
      x = c(xval, rev(xval)),
      y = c(summary, rev(recent)),
      col = colors[[c]]
    )
    recent <- summary
  }
}

do_model1 <- function(dom, norm_perct = NULL, ma = 0) {
  if (length(dim(dom)) != 5)
    stop("Dom dimensions must be equal to 5 [lon, lat, level, time]")
  s <- colSums(dom, dims = 2)
  w <- uv2hd(s)
  if (!is.null(norm_perct)) # Normalize each level seperately
    w[,, "h"] <- apply(w[,, "h"], 1, normalize, perct = norm_perct)
  if (ma != 0) {
    ma_d <- t(apply(w[,, "d"], 1, SMA, n = ma))
    ma_h <- t(apply(w[,, "h"], 1, SMA, n = ma))
    dimnames(ma_d) <- dimnames(w[,, "d"])
    dimnames(ma_h) <- dimnames(w[,, "h"])
    w[,, "d"] <- ma_d
    w[,, "h"] <- ma_h
  }
  d_comp  <- deg_to_comp(w[,, "d"])
  d_comp2 <- deg_to_comp(w[,, "d"], bin_names = seq(1, 16))
  # w <- round(w, 3)
  return(list(w = w, d_comp = d_comp, d_comp2 = d_comp2))
}

do_model2 <- function(dom, alfa = 0, ma=0, norm_perct = NULL) {
  s <- colSums(dom, dims = 2)
  rs <- rotax(s, alfa = alfa)
  if (ma != 0) {
    ma_u <- t(apply(rs[,, 1], 1, SMA, n = ma))
    ma_v <- t(apply(rs[,, 2], 1, SMA, n = ma))
    dimnames(ma_u) <- dimnames(rs)[1:2]
    dimnames(ma_v) <- dimnames(rs)[1:2]
   rs[,, 1] <- ma_u
   rs[,, 2] <- ma_v
  }
  if (!is.null(norm_perct))
    rs <- aperm(apply(rs, c(1, 3), normalize, perct = norm_perct), c(2, 1, 3))
  return( round(rs, 3) )
}

check_domain <- function(x, pm, rotate_degree = 0) {
  # focus on 1000 hPa
  x <- x[,,1,,]
  # pm <- as.vector(read_pm10()[, 4])
  fl_cor_alfa <- file.path(dir_data, "cor_alfa.rdata")
  if (file.exists(fl_cor_alfa)) {
    load(fl_cor_alfa)
  } else {
    cor_alfa <- calc_cor(x, pm, alfa = seq(0, 360, 1))
    save(cor_alfa, file = fl_cor_alfa)
  }
  max_ind <-  which(cor_alfa == max(cor_alfa), arr.ind = T)[1,]
  print(max_ind)
  # MAX. Correlation = 0.460589
  # On level = 1000 hPa at lon = 30 and lat = 42.5
  # At alfa = 9 degree for v component
  # At alfa = 279 degree for u component
  comp <- "u"
  lon <- "30"; lat <- "42.5"
  u   <- cor_alfa[, 1, lon, lat ]
  v   <- cor_alfa[, 2, lon, lat ]
  x   <- as.numeric(names(u))
  plot(x, u, type = "l", las = 1, col = "green", lty = 1,
       xlab = "alfa", ylab = "Correlation",
       main = "cor(pm10 ~ u) & cor(pm10 ~ v) \n variability by alfa",
       sub = paste0("lon=", lon, " lat=", lat))
  lines(x, v, col = "red")
  abline(h = 0, col = "blue", lty = 2)

  min_ind <- which(cor_alfa == min(cor_alfa), arr.ind = T)[1,]
  print(min_ind)
  # MIN. Correlation = -0.460589
  # On level = 1000 hPa at lon = 30 and lat = 42.5
  # At alfa = 189 degree for v component
  # At alfa = 99 degree for u component

  # I'm interested in 1000 hPa level and
  # the grid point lon = 30 and lat = 42.5
  # at alfa = 9 degree.
  alfa <- "9"; comp <- "v"
  lon <- "30"; lat <- "42.5"
  x  <- cor_alfa[alfa, comp, lon, lat ]
  y  <- as.numeric(names(x))
  y1 <- 1013.6 - y
  plot(x, y1, col = "red", type = "l", yaxt = "n",
       xlab = "Correlation", ylab = "Levels",
       main = "pm10 ~ v(alfa=9) & pm10 ~ u(alfa=299) \n variability by level",
       sub = paste0("lon=", lon, " lat=", lat))
  y_lab <- as.character(y)
  y_lab[13:16] <- ""
  axis(2, at = y1, labels = y_lab, las = 2, cex = 0.5)
  points(x, y1, col = "red", pch = 16)

  alfa <- "299"
  x  <- cor_alfa[alfa, "u", lon, lat, ]
  y  <- as.numeric(names(x))
  y1 <- 1013.6 - y
  lines(x, y1, col = "blue")
  points(x, y1, col = "blue", pch = 16)

  grid()

  cor_v_alfa <- cor_alfa[, "v", lon, lat, "1000"]
  cor_u_alfa <- cor_alfa[, "u", lon, lat, "1000"]
  degrees <- as.numeric(names(cor_v_alfa))
  radial.lim <- c(-0.5, 0.5)
  polar.plot(cor_v_alfa, degrees,
             line.col = "red", show.grid.labels = 6, poly.col = "pink",
             point.symbols = 16, point.col = "red", radial.lim = radial.lim,
             rp.type = "p", start = 90, clockwise = T)
  polar.plot(cor_v_alfa[10], degrees[10], add = T,
             line.col = "red", show.grid.labels = 6, poly.col = "pink",
             point.symbols = 16, point.col = "red", radial.lim = radial.lim,
             rp.type = "s", start = 90, clockwise = T, radlab = T,
             labels = NA, label.pos = degrees[10])
  polar.plot(cor_u_alfa, degrees,
             line.col = "blue", poly.col = "cyan",
             point.symbols = 16, point.col = "blue", radial.lim = radial.lim,
             rp.type = "p", start = 90, clockwise = T, add = T)
  polar.plot(cor_u_alfa[300], degrees[300], add = T,
             line.col = "blue", show.grid.labels = 6, poly.col = "cyan",
             point.symbols = 16, point.col = "blue", radial.lim = radial.lim,
             rp.type = "s", start = 90, clockwise = T)

  alfa <- "0"; comp <- "v"
  lat <- 21:1
  for (l in 1:17) {
    level <- dimnames(cor_alfa)[[5]][l]
    mat <- cor_alfa[alfa, comp, ,lat,l]
    filled.contour(x = as.numeric(rownames(mat)),
                   y = as.numeric(colnames(mat)),
                   mat,
                   main = paste0("Level= ", level,
                                 " hPa (pm10 ~ Wind-", comp, " component)"),
                   sub = paste0("alfa=", alfa),
                   plot.axes = {
                     axis(1); axis(2);
                     map(add = TRUE)
                   })
  }
}

lineer_model <- function(formula, col = "blue", bg = "cornflowerblue",
                         angle = 45, zlab = NULL, ...) {
  mdl <- lm(formula)
  if (ncol(mdl$model) == 3) {
    if (require(scatterplot3d)) {
      scatterplot3d(mdl$model[3:1], pch = 21, angle = angle, color = "blue",
                    zlab = zlab,
                    bg = "cornflowerblue", ...)$plane3d(mdl, lty.box = "solid",
                                                        col = "red")
      return(mdl)
    }
  }
  plot(mdl$model[, 2:1], pch = 21, col = col, bg = bg, ...)
  abline(mdl, col = "red", lwd = 2)
  return(mdl)
}


moving_avg <- function(x, ma = seq(0, 14, 1), dim = 2L) {
  dm <- dim(x); dn <- dimnames(x)
  x2 <- na.locf(as.matrix(x), fromLast = T)
  mavs <- vapply(ma,
                 function(n) if (n != 0) apply(x2, dim, SMA, n = n) else x,
                 outer(1:dm[1], 1:dm[2]))
  dim(mavs) <- c(dm, ma = length(ma))
  dimnames(mavs) <- append(dn, list(ma = ma))
  return(mavs)
}

acf_test <- function() {
  acfs <- vapply(1:dim(mavs)[3],
                 function(i) {
                   a <- acf(mavs[, c(1,5), i], lag.max = 20,
                            na.action = na.pass, plot = F)
                   a_max <- t(sapply(1:length(a$snames), function(j) cbind(which(a$acf[, -j, j] == max(a$acf[, -j, j]), arr.ind = T), max(a$acf[, -j, j])) )) + c(-1,1)
                   return(a_max)
                 }, outer(1:5,1:3))

  l <- list()
  for (i in 1:dim(mavs)[3]) {
    mav <- mavs[,,i]
    for (j in 1:(dim(mav)[2] - 1)) {
      for (k in (j + 1):dim(mav)[2]) {
        a <- ccf(mav[, j], mav[, k], lag.max = 20, na.action = na.pass, plot = F)
        acf_max <- max(a$acf)
        acf_max_lag <- a$lag[a$acf == acf_max]
        l <- append(l, list(data.frame(i, colnames(mav[,j, drop = F]), colnames(mav[,k, drop = F]), acf_max, acf_max_lag)))
      }
    }
  }

  # Check correlations by moving average and lag
  # this shows yo moving average anf applying lag has a meaning or not.
  ma_steps <- seq(0, 14, 1)
  mavs <- array(0, dim = c(time = nrow(mwp), vars = ncol(mwp),
                           ma = length(ma_steps)))
  dimnames(mavs) <- append(dimnames(mwp), list(ma = ma_steps))
  result <- array(0, dim = c(ma = length(ma_steps), measurements = 3))
  dimnames(result) <- list(ma = ma_steps,
                           measurement = list("cor", "acf_max", "acf_lag"))
  for (ma in ma_steps) {
    mav <- if (ma != 0) apply(mwp, 2, SMA, n = ma) else mwp
    cr  <- cor(mav[, 2], mav[, 5], use = "pairwise.complete.obs")
    a   <- ccf(mav[, 2], mav[, 5], lag.max = 100, na.action = na.pass, plot = F)
    acf_max <- max(a$acf)
    acf_max_lag <- a$lag[a$acf == acf_max]
    # acf_max <- round(acf_max, 3)
    result[ma + 1, ] <- c(cr, acf_max, acf_max_lag)
    mavs[,, ma + 1] <- mav
  }
  cat("max cor = ", max(result[, 1]), " and on ma = ",
      which(result[, 1] == max(result[, 1])) - 1, "\n")
  print(result)
}

check_relations <- function(pm, hgt, m) {
  x <- cbind(as.matrix(pm), as.matrix(hgt), as.matrix(m[1,,]))
  vars <- colnames(x)
  names(dim(x)) <- c("time", "vars")
  dimnames(x) <- list(time = colnames(m), vars = vars)

  cor(x, use = "pairwise.complete.obs")

  mavs <- moving_avg(x) # calc moving averages
  crs <- vapply(1:dim(mavs)[3],
                function(i) cor(mavs[,,i], use = "pairwise.complete.obs"),
                outer(1:dim(mavs)[2], 1:dim(mavs)[2]))
  dn <- append(dimnames(crs)[1:2], list(ma = dimnames(mavs)[[3]]))
  dim(crs) <- c(dim(mavs)[2], dim(mavs)[2:3])
  dimnames(crs) <- dn



  # regression
  lm_pm     <- as.vector(mavs[, 1, 4])
  lm_pm_log <- as.vector(mavs[, 2, 4])
  lm_hgt1   <- as.vector(mavs[, 3, 4])
  lm_u      <- as.vector(mavs[, 4, 4])
  lm_v      <- as.vector(mavs[, 5, 4])

  lm_pm     <- as.vector(x[, 1])
  lm_pm_log <- as.vector(x[, 2])
  lm_hgt1   <- as.vector(x[, 3])
  lm_u      <- as.vector(x[, 4])
  lm_v      <- as.vector(x[, 5])

  mdl <- lineer_model(lm_pm_log ~ lm_hgt1,
                      type = "p", xlab = "hgt1",
                      ylab = "log(pm10)")

  mdl <- lineer_model(lm_pm_log ~ lm_v,
                      type = "p", xlab = "Wind-v",
                      ylab = "log(pm10)")



  fit <- reg_tree(x[,-1], main = "(pm10 ~ HGHT1 + u + v) for alfa = 9")
  plot_reg_tree(fit, 3)
  plot_reg_tree(fit, 4)
  plot_reg_tree(fit, 6)
  plot_reg_tree(fit, 7)

  df <- x[, c(1, 5)]
  rownames(df) <- 1:nrow(df)
  df[, 1] <- normalize(df[, 1])
  df[, 2] <- normalize(df[, 2])

  df <- df[0:50, ]
  df2 <- data.frame(value = as.vector(df), index = as.numeric(rownames(df)),
                    series = rep(colnames(df), each = nrow(df)))
  ggplot(df2, aes(x = index, y = value, fill = series)) +
         geom_area(position = "stack")

}

reg_tree <- function(x, ...) {
  d <- x
  if (!is.formula(d)) {
    stopifnot(!is.null(dim(d)) & length(dim(d)) == 2)
    d <- as.data.frame(d)
  }
  fit <- rpart(d, method = "anova")
  # printcp(fit)
  # plotcp(fit)
  # summary(fit)
  par(xpd = TRUE)
  plot(fit, uniform = TRUE, col = "gray", ...)
  text(fit, use.n = T, all = T, pos = 1,
       fancy = T, fheight = 4, fwidth = 7, bg = "antiquewhite")
  return(fit)
}

plot_reg_tree <- function(fit, split_id) {
  dates <- names(fit$where[fit$where == split_id[1]])
  for (i in split_id[-1]) {
    dates <- c(dates, names(fit$where[fit$where == i]))
  }
  # dates <- c("2008-04-03" , "2008-11-03", "2009-01-29", "2009-11-25", "2011-02-03", "2011-11-29")
  # boxplot(as.matrix(pm[dates,]))
  # head(pm)
  # filter ua and va by dates
  ua <- w[,,1,,1]
  va <- w[,,1,,2]
  u <- aperm(ua[,,dates], c(3,1,2))
  v <- aperm(va[,,dates], c(3,1,2))
  u_sd <-apply(u, c(2,3), sd, na.rm = T)
  v_sd <-apply(v, c(2,3), sd, na.rm = T)
  u_med <- apply(u, c(2,3), median, na.rm = T)
  v_med <- apply(v, c(2,3), median, na.rm = T)
  u_max <- apply(u, c(2,3), max, na.rm = T)
  v_max <- apply(v, c(2,3), max, na.rm = T)

  u <- colMeans(u); v <- colMeans(v)

  # plot_streamlines(u - u_sd, v - v_sd, main = paste0(as.character(split_id), " -SD", collapse = " "))
  plot_streamlines(u, v, main = paste0(as.character(split_id), " Mean", collapse = " "))
  # plot_streamlines(u_med, v_med, main = paste0(as.character(split_id), " Median", collapse = " "))
  # plot_streamlines(u_max, v_max, main = paste0(as.character(split_id), " Max", collapse = " "))
  # plot_streamlines(u + u_sd, v + v_sd, main = paste0(as.character(split_id), " +SD", collapse = " "))
  # par(xpd = T)
  # plot_vector_field(u, v, main = paste0(as.character(split_id), collapse = " "))
}

plot_streamline_day <- function(dates) {
  lp <- list()
  for (d in dates) {
    u <- ua[,,1,d]
    v <- va[,,1,d]
    p <- plot_streamlines(u, v, main = d)
    lp <- append(lp, list(p))
    print(p)
  }
  return(lp)
}

plot.stack <- function(x, y) {

  # plot stacked line
  df <- mwp[, c(1, 2)]
  rownames(df) <- 1:nrow(df)
  df[, 1] <- normalize(df[, 1])
  df[, 2] <- normalize(df[, 2])

  df <- df[750:900, ]
  df2 <- data.frame(value = as.vector(df), index = as.numeric(rownames(df)),
                    series = rep(colnames(df), each = nrow(df)))
  ggplot(df2, aes(x = index, y = value, fill = series)) +
         geom_area(position = "stack")
}

create_domain <- function(w, lon = 30.0, lat = 42.5) {
  lons <- as.character(lon)
  lats <- as.character(lat)
  ilons <- as.numeric(factor(lons, levels = dimnames(w)[[1]]))
  ilats <- as.numeric(factor(lats, levels = dimnames(w)[[2]]))

  selection_array <- w * 0
  for (i in 1:length(ilons))
    selection_array[ilons[i], ilats[i],,, ] <- 1

  w <- w * selection_array
  # remove zero rows
  remove_zeros <- function(x) {
    x <- x[rowSums(abs(x)) != 0,,,,, drop = F]
    x <- aperm(x, c(2, 1, 3, 4, 5))
    x <- x[rowSums(abs(x)) != 0,,,,, drop = F]
    return(aperm(x, c(2, 1, 3, 4, 5)))
  }
  w <- remove_zeros(w)
  dn <- dimnames(w)
  names(dim(w)) <- names(dn)
  dimnames(w) <- dn
  return(w)
}

load_wind <- function() {
  fl <- "wind.rdata"
  pf <- file.path(dir_data, fl)
  if (!file.exists(pf)) {
    u_nc <- nc_open(file.path(dir_data, "uwnd.nc"))
    v_nc <- nc_open(file.path(dir_data, "vwnd.nc"))
    u    <- get_data(u_nc)
    v    <- get_data(v_nc)
    nc_close(u_nc)
    nc_close(v_nc)
    w <- array(c(u, v), dim = c(dim(u), component = 2))
    dimnames(w) <- append(dimnames(u), list(component = c("u", "v")))
    save_rdata(w, file = fl)
    cat("Loaded *wnd.nc\n")
  }else{
    load(pf)
    cat("Loaded", fl, "\n")
  }
  return(w)
}

load_interpolated_wind_1000 <- function(x) {
  nv <- deparse(substitute(x))
  fl <- "interpolated_1000.rdata"
  pf <- file.path(dir_data, fl)
  if (!file.exists(pf)) {
    w1 <- x[,,1,,]
    dim(w1) <- dim(x)[c(1, 2, 4, 5)]
    dimnames(w1) <- dimnames(x)[c(1, 2, 4, 5)]
    cat("Calculationg interpolation for", nv, "\n")
    system.time( wi <- interp2(w1, dx = 0.2) )
    save_rdata(wi, file = fl)
  } else {
    load(pf)
    cat("Loaded interpolated wind for", nv, "\n")
  }
  return(wi)
}

get_correlations <- function(x, pm) {
  nv1 <- deparse(substitute(x))
  nv2 <- deparse(substitute(pm))
  fl <- paste0("cor_alfa_", nv1, "_", nv2, ".rdata")
  pf <- file.path(dir_data, fl)
  if (!file.exists(pf)) {
    cat("Calculating correlations for", nv1, "~", nv2, "\n")
    t <- system.time( cor_alfa <- calc_cor(x, pm, alfa = seq(0, 360, 1)) )
    print(t)
    save_rdata(cor_alfa, file = fl)
  } else {
    load(pf)
    cat("Loaded correlations for", nv1, "~", nv2, "\n")
  }
  return(cor_alfa)
}
