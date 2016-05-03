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
library(lubridate)
library(PCICt)
library(maps)
library(plotrix)
library(ncdf4)
library(ncdf4.helpers)
library(TTR)
library(timeSeries)
library(parallel)
library(akima)
library(rasterVis)

# INITIALS
dir_data <- "data"
dir_out  <- "output"

# download map
dir_countries <- "countries_shp"
fl <- "ne_50m_admin_0_countries"
if (!dir.exists(dir_countries)) {
  fle <- paste0(fl, ".zip")
  url <- paste0("http://naciscdn.org/naturalearth/50m/cultural/", fle)
  download.file(url, fle, method = "curl")
  unzip(fle, exdir = dir_countries)
}
cntry <- readOGR(dsn = paste0("./", dir_countries),
                 layer = "ne_50m_admin_0_countries", stringsAsFactors = TRUE)
rm(fl, dir_countries)

# http://stackoverflow.com/questions/14500707/select-along-one-of-n-dimensions-in-array
index_array <- function(x, dim, value, drop = T) {
  # Create list representing arguments supplied to [
  # bquote() creates an object corresponding to a missing argument
  indices <- rep(list(bquote()), length(dim(x)))
  indices[[dim]] <- value

  # Generate the call to [
  call <- as.call(c(
    list(as.name("["), quote(x)),
    indices,
    list(drop = drop)))
  # Print it, just to make it easier to see what's going on
  # print(call)

  # Finally, evaluate it
  eval(call)
}

save_rdata <- function(..., file = stop("'file' must be specified")) {
  if (nchar(file_ext(file)) == 0) file <- paste0(file, ".rdata")
  save(..., file = file.path(dir_data, file), envir = parent.frame())
}

savecsv <- function(x, file = "") {
  dts <- if (is.xts(x)) index(x) else as.POSIXct(rownames(x), tz = "GMT")
  dts <- format(dts, "%d.%m.%Y")
  df <- as.data.frame(x)
  df <- cbind(dts, df)
  colnames(df)[[1]] <- "date"
  write.table(df, file = file, quote = F, sep = ";", row.names = F)
}

read_pm10 <- function(f = "NightResults12042016.csv", shift = 0) {
  d <- read.table(file.path(dir_data, f))
  # find PM
  indices <- c(5, which(grepl("PM", colnames(d)) == T))
  d <- d[, indices] # filter only PM data
  d <- d[, colSums(is.na(d)) < nrow(d)] # remove full NA cols
  d[, 1] <- as.POSIXct(d[, 1], tz = "GMT") + days(shift)
  d <- xts(d[, -1], order.by = d[, 1])

  dts <- seq(as.POSIXct("2008-01-01", tz = "GMT"),
             as.POSIXct("2015-12-27", tz = "GMT"), by = "days")
  # dts[which(dts %in% index(d) == F)]
  d <- merge(d, zoo(, dts), all = T)
  d <- d["2008/2012"]
  return(d)
}

read_wind <- function(f="NightResults12042016.csv", shift=0) {
  d <- read.table(file.path(dir_data, f))
  # find wind components
  indices <- c(5, which(grepl("Wind", colnames(d)) == T))
  d <- d[, indices] # filter only wind data
  d <- d[, colSums(is.na(d)) < nrow(d)] # remove full NA cols
  d[, 1] <- as.POSIXct(d[, 1], tz = "GMT") + days(shift)
  d <- xts(d[, -1], order.by = d[, 1])

  dts <- seq(as.POSIXct("2008-01-01", tz = "GMT"),
             as.POSIXct("2015-12-27", tz = "GMT"), by = "days")
  # dts[which(dts %in% index(d) == F)]
  d <- merge(d, zoo(, dts), all = T)
  d <- d["2008/2012"]
  dir <- d[, seq(1, ncol(d), 2)]
  hs <- d[, seq(2, ncol(d), 2)]
  return(list(h = hs, d = dir))
}

read_hgt <- function(f="NightResults12042016.csv", shift=0) {
  d <- read.table(file.path(dir_data, f))
  # fin wind components
  indices <- c(5, which(grepl("HGHT", colnames(d)) == T))
  d <- d[, indices] # filter only wind data
  d <- d[, colSums(is.na(d)) < nrow(d)] # remove full NA cols
  d[, 1] <- as.POSIXct(d[, 1], tz = "GMT") + days(shift)
  d <- xts(d[, -1], order.by = d[, 1])

  dts <- seq(as.POSIXct("2008-01-01", tz = "GMT"),
             as.POSIXct("2015-12-27", tz = "GMT"), by = "days")
  # dts[which(dts %in% index(d) == F)]
  d <- merge(d, zoo(, dts), all = T)
  d <- d["2008/2012"]
  return(d)
}

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


deg_to_comp <- function(x, bin_names=c("N", "NNE", "NE", "ENE",
                                       "E", "ESE", "SE", "SSE",
                                       "S", "SSW", "SW", "WSW",
                                       "W", "WNW", "NW", "NNW")) {
  x[x == 0] <- 360
  i <- 360 / length(bin_names)
  val <- round( (x / i) + .5)
  return(matrix(bin_names[val], nrow(x), ncol(x), dimnames = dimnames(x)))
}


uv_to_hd <- function(u, v) {
  if (all(length(u) != length(v))) stop("u and v length MUST be the same.")
  h <- sqrt(u ^ 2 + v ^ 2)
  d <- ((atan2(u / h, v / h) * 180 / pi) + 180) %% 360
  return(list(h = h, d = d))
}


hd_to_uv <- function(h, d) {
  if (all(length(h) != length(d))) stop("'h' and 'd' length MUST be the same.")
  h <- h * 0.514444444 # convert knots to m/s
  d <- (d - 180) * pi / 180
  return(list(u = h * sin(d), v = h * cos(d)))
}

# Convert u,v wind components to horizontal speed and direction
uv2hd <- function(..., drop = T, degree = T, reverse = T) {
  l <- list(...)
  if (length(l) > 2) {
    x <- do.call(c, l)
  }else if (length(l) == 2) {
    if (all(dim(l[[1]]) == dim(l[[2]]))) {
      dm <- dim(l[[1]])
      x <- array(c(l[[1]], l[[2]]), dim = c(dm, 2))
    }
  } else {
    x <- if (is.array(l[[1]])) l[[1]] else do.call(c, l)
  }
  if (is.vector(x) && (length(x) %% 2) == 0)
    x <- array(x, dim = c(length(x) / 2, 2))

  comp_loc <- which(dim(x) == 2) # find uv dim order
  if (length(comp_loc) == 0)
    stop("x does not have a dim = 2 represents u and v")
  dm <- dim(x)
  dmn <- dimnames(x)
  l <- length(dim(x))
  r <- aperm(x, c(comp_loc, (1:l)[-comp_loc])) # get uv dim to start
  h <- sqrt(colSums(r ^ 2))
  u <- index_array(r, 1, 1)
  v <- index_array(r, 1, 2)
  num_pi <- if (degree) 180 else pi
  add_pi <- if (reverse) num_pi else 0
  d <- ((atan2(u / h, v / h) * num_pi / pi) + add_pi)
  if (degree) d <- d  %% 360

  r <- array(c(h, d), dim = dm)
  dmn[[comp_loc]] <- c("h", "d")
  dimnames(r) <- dmn

  r[is.nan(r)] <- 0
  if (drop) r <- drop(r)
  return(r)
}

hd2uv <- function(..., drop = T, degree = T) {
  l <- list(...)
  if (length(l) >= 2) {
    x <- do.call(c, l)
  } else {
    x <- if (is.array(l[[1]])) l[[1]] else do.call(c, l)
  }

  if (is.vector(x) && (length(x) %% 2) == 0)
    x <- array(x, dim = c(length(x) / 2, 2))

  comp_loc <- which(dim(x) == 2) # find hd dim order
  if (length(comp_loc) == 0)
    stop("x does not have a dim = 2 represents h and d")
  dm <- dim(x)
  dmn <- dimnames(x)
  l <- length(dim(x))
  r <- aperm(x, c(comp_loc, (1:l)[-comp_loc])) # get uv dim to start
  h <- index_array(r, 1, 1)
  d <- index_array(r, 1, 2)
  # d <- (d - 180) * pi/180
  num_pi <- if (degree) 180 else pi
  d <- (d - num_pi) / num_pi
  u <- h * sinpi(d)
  v <- h * cospi(d)
  r <- array(c(u, v), dim = dm)
  dmn[[comp_loc]] <- c("u", "v")
  dimnames(r) <- dmn
  r[is.nan(r)] <- 0
  if (drop) r <- drop(r)
}

get_data <- function(nc, lats=NULL, lons=NULL, levels=NULL, times=NULL) {
  c_to_i <- function(dim.name, vals) {
    return(match(vals, ncvar_get(nc, dim.name)))
  }
  i_to_c <- function(dimname, indices) {
    return(ncvar_get(nc, dimname)[indices])
  }

  var_name <- nc$var[[2]]$name

  ilats <- if (is.null(lats)) 1:nc$dim$lat$len else c_to_i("lat", lats)
  lats  <- i_to_c("lat", ilats)

  ilons <- if (is.null(lons)) 1:nc$dim$lon$len else c_to_i("lon", lons)
  lons  <- i_to_c("lon", ilons)

  ts <- nc.get.time.series(nc)
  if (is.character(times)) itimes <- which(as.character(ts) == times)
  if (is.null(times)) itimes <- 1:nc$dim$time$len
  times <- as.character(ts[itimes])

  if (is.null(levels)) ilevels <- 1:nc$dim$level$len
  if (is.character(levels)) ilevels <- c_to_i("level", levels)
  levels <- i_to_c("level", ilevels)

  data <- nc.get.var.subset.by.axes(nc, var_name,
                                   list(lon   = ilons,
                                        lat   = ilats,
                                        level = ilevels,
                                        time  = itimes),
                                   nc.get.dim.names(nc)[1:4])
  dimnames(data) <- list(lon = lons, lat = lats, level = levels, time = times)
  return(data)
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

plot_vector_field <- function(u, v, lon=NULL, lat=NULL, show_grids=T,
                     add=F, main="") {
  if (is.null(lon)) lon <- as.numeric(rownames(u))
  if (is.null(lat)) lat <- as.numeric(colnames(u))
  g <- create_grids(lon, lat)

  if (!add) {
    par(mar = c(0, 0, 0, 0))
    map("world", xlim = range(g[, 1]), ylim = range(g[, 2]),
        col = "gray95", fill = T)
    if (show_grids) points(g, pch = 19, col = "red", cex = 0.2)
  }
  w <- uv_to_hd(u, v)
  vectorField( (270 - w$d), w$h, g[, 1], g[, 2], scale = 1.8,
               vecspec = "deg", col = "gray50")
  if (main != "") mtext(main, 3, line = -6)
}

plot_vector <- function(ua, va, dom=NULL, model=NULL, level=1, time=1,
                        show_vectors=T, show_grids=T) {

  ilev  <- if (is.character(level)) which(dimnames(ua)[[3]] == level) else level
  level <- dimnames(ua)[[3]][ilev]

  itime <- if (is.character(time)) which(dimnames(ua)[[4]] == time) else time
  time  <- dimnames(ua)[[4]][itime]

  u <- ua[,, ilev, itime]
  v <- va[,, ilev, itime]
  plot_vector_field(u, v, main = paste0(time, " / ", level, " hPa"))

  if ( !is.null(dom) ) {
    grids <- cbind(as.numeric(dom$lons), as.numeric(dom$lats))
    hpts  <- chull(grids)
    hpts  <- c(hpts, hpts[1])
    polygon(grids[hpts, ], col = rgb(0, 0, 1, 0.2))
    points(grids, pch = 19, col = "blue", cex = 0.3)
    #
    if ( !is.null(model) ) {
      vectorField( (270 - model$w$d[itime, ilev]), model$w$h[itime, ilev],
                   xpos = c(22.5, 25), ypos = c(35, 37.5), vecspec = "deg",
                   scale = 2, col = "red" )
    }
  }
}


interpolate <- function(u, v) {
  u <- u[,ncol(u):1]
  v <- v[,ncol(v):1]
  x <- as.numeric(dimnames(u)$lon)
  y <- as.numeric(dimnames(u)$lat)
  ui <- bicubic.grid(x, y, u, xlim = range(x), ylim = range(y), dx = 0.1, dy = 0.1)
  vi <- bicubic.grid(x, y, v, xlim = range(x), ylim = range(y), dx = 0.1, dy = 0.1)

  u <- ui$z[, ncol(ui$z):1]
  dimnames(u) <- list(lon = ui$x, lat = ui$y[length(ui$y):1])
  v <- vi$z[, ncol(vi$z):1]
  dimnames(v) <- list(lon = vi$x, lat = vi$y[length(vi$y):1])
  return(list(u = u, v = v))
}

calc_interp <- function(m, dx = NULL, dy = dx) {
  stopifnot( length(dim(m)) == 2 )
  rev_x <- rev_y <- F

  x <- as.numeric(rownames(m))
  y <- as.numeric(colnames(m))
  if (length(x) == 0) x <- 1:nrow(m)
  if (length(y) == 0) y <- 1:ncol(m)

  if (is.null(dx)) dx <- abs((x[2] - x[1]) / 2)
  if (is.null(dy)) dy <- abs((y[2] - y[1]) / 2)
  new_x <- seq(min(x), max(x), dx)
  new_y <- seq(min(y), max(y), dy)
  new_nr <- length(new_x)
  new_nc <- length(new_y)

  if (!all(x == cummax(x))) {
    x <- x[length(x):1]
    m <- m[nrow(m):1,]
    new_x <- rev(new_x)
    rev_x = T
  }
  if (!all(y == cummax(y))) {
    y <- y[length(y):1]
    m <- m[,ncol(m):1]
    new_y <- rev(new_y)
    rev_y = T
  }

  return(list(x = x, y = y, m = m,
              dx = dx, dy = dy,
              new_x = new_x, new_y = new_y,
              new_nr = new_nr, new_nc = new_nc,
              rev_x = rev_x, rev_y = rev_y))
}

interpb <- function(m, dx = NULL, dy = dx) {
  stopifnot( length(dim(m)) == 2 )
  r <- calc_interp(m, dx, dy)
  rm(m)
  library(akima)
  mi <- bicubic.grid(r$x, r$y, r$m,
                     xlim = range(r$x), ylim = range(r$y),
                     dx = r$dx, dy = r$dy)

  if (r$rev_x) {
    mi$x <- rev(mi$x)
    mi$z <- mi$z[nrow(mi$z):1,]
  }
  if (r$rev_y) {
    mi$y <- rev(mi$y)
    mi$z <- mi$z[,ncol(mi$z):1]
  }
  rm(r)
  names(dim(mi$z)) <- c("lon", "lat")
  dimnames(mi$z) <- list(lon = mi$x, lat = mi$y)
  return(mi$z)
}

interp2 <- function(x, dx = NULL, dy = dx) {
  stopifnot( length(dim(x)) > 1 )
  dm <- dim(x)[1:2]
  a <- array(as.vector(x)[1:(dm[1] * dm[2])], dim = dm)
  dim(a) <- dm; dimnames(a) <- dimnames(x)[1:2]
  r <- calc_interp(a, dx, dy)
  dm <- (1:length(dim(x)))[-(1:2)]
  if (length(dm) > 0) {
    cl <- makeCluster(detectCores())
    clusterExport(cl, c("interpb", "calc_interp", "r"), envir = environment())
    mi <- parApply(cl, x, dm, function(m) interpb(m, r$dx, r$dy))
    stopCluster(cl)
    mi <- array( mi, dim = c(r$new_nr, r$new_nc, dim(x)[-(1:2)]) )
    names(dim(mi)) <- names(dim(x))
    dimnames(mi) <- append(list(lon = r$new_x, lat = r$new_y), dimnames(x)[-(1:2)])
  } else {
    mi <- interpb(x, r$dx, r$dy)
  }
  return(mi)
}

plot_vector_field2 <- function(u, v, ...) {
  # proj <- CRS('+proj=longlat +datum=WGS84')
  l <- interpolate( u, v)
  u <- l$u; v <- l$v
  lon <- as.numeric(rownames(u))
  lat <- as.numeric(colnames(u))
  xmn = min(lon); xmx = max(lon)
  ymn = min(lat); ymx = max(lat)
  u <- t(u); v <- t(v)
  u1 <- raster(x = u, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx)
  v1 <- raster(x = v, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx)
  w <- brick(u1, v1)
  projection(w) <- CRS("+init=epsg:4326")
  extent(w) <- c(xmn, xmx, ymn, ymx)
  slope <- sqrt(w[[1]] ^ 2 + w[[2]] ^ 2)
  print(max(slope@data@values))
  # aspect <- atan2(w[[1]], w[[2]])
  vectorplot(w * 2, isField = "dXY", region = slope, margin = F, lwd = 1.2,
             par.settings = RdBuTheme(), narrows = 1000, at = 0:10, ...) +
    layer(sp.polygons(cntry, fill = 'transparent', col = "blue", alpha = 0.5))
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

rotax <- function(..., alfa=45, degree=T, right=T) {
  l <- list(...)
  if (length(l) > 2) {
    x <- do.call(c, l)
  }else if (length(l) == 2) {
    if (all(dim(l[[1]]) == dim(l[[2]]))) {
      dm <- dim(l[[1]])
      x <- array(c(l[[1]], l[[2]]), dim = c(dm, 2))
    }
  } else {
    x <- if (is.array(l[[1]])) l[[1]] else do.call(c, l)
  }
  if (is.vector(x) && (length(x) %% 2) == 0)
    x <- array(x, dim = c(length(x) / 2, 2))
  comp_loc <- which(dim(x) == 2) # find uv dim order
  if (length(comp_loc) == 0)
    stop("x does not have a dim = 2 represents u and v")
  dm <- dim(x); dn <- dimnames(x)
  u <- index_array(x, comp_loc, 1)
  v <- index_array(x, comp_loc, 2)
  a <- if (degree) alfa / 180 else alfa / pi
  c <- cos(pi*a); s <- sin(pi*a)
  if (right) {
    un <- u * c - v * s
    vn <- u * s + v * c
  } else {
    un <- u * c + v * s
    vn <- v * c - u * s
  }
  x <- array(c(un, vn), dim = dm)
  dimnames(x) <- dn
  return(drop(x))
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

calc_cor <- function(x, pm, alfa = seq(0, 360, 20), par = T) {
  # pm <- as.vector(pm)
  calc <- function(uv) {
    rot <- function(a) rotax(uv, alfa = a)
    cor2 <- function(x, y) cor(y, x, use = "pairwise.complete.obs")
    r <- vapply(alfa, rot, outer(1:dim(uv)[1], 1:dim(uv)[2]))
    ret <- t(apply(r, c(2, 3), cor, y = pm, use = "pairwise.complete.obs"))
    # print(max(ret))
    return(ret)
  }
  comp_loc <- which(dim(x) == 2) # find uv dim order
  if (length(comp_loc) == 0)
    stop("x does not have a dim = 2 represents u and v")
  time_loc <- which(dim(x) == nrow(pm)) # find time dim order
  if (length(time_loc) == 0)
    stop("x does not have a dim suitable with pm")
  if (par) {
    # cl <- makeCluster(4, type = "FORK")
    cl <- makeCluster(detectCores() - 1, outfile = "")
    clusterExport(cl, c("pm", "index_array", "rotax", "alfa"), envir = environment())
    r <- parApply(cl, x, (1:length(dim(x)))[-c(comp_loc, time_loc)], calc)
    stopCluster(cl)
  } else {
    r <- apply(x, (1:length(dim(x)))[-c(comp_loc, time_loc)], calc)
  }
  r <- array(r, dim = c(length(alfa), 2, dim(x)[-c(comp_loc, time_loc)]))
  names(dim(r)) <- c("alfa", "component", names(dim(x)[-c(comp_loc, time_loc)]))
  dimnames(r) <- append(list(alfa = alfa, component = c("u", "v")), dimnames(x)[-c(comp_loc, time_loc)])
  dim_order <- 1:length(dim(r))
  comp_loc <- which(dim(r) == 2)
  alfa_loc <- which(dim(r) == length(alfa))
  dim_comp <- dim_order[comp_loc]
  dim_alfa <- dim_order[alfa_loc]
  dim_remain <- dim_order[-c(comp_loc, alfa_loc)]
  r <- aperm(r, c(dim_remain, dim_comp, dim_alfa))
  return(r)
}

stat_cor <- function(x) {
  get_names <- function(dn, i) {
    m <- vector()
    for (j in 1:length(i)) m <- c(m, dn[[names(i[j])]][i[j]])
    return(m)
  }
  dn <- dimnames(x)
  max_cor <- max(x); min_cor <- min(x)
  imax <-  which(x == max_cor, arr.ind = T)[1,]
  max <- get_names(dn, imax)
  imin <- which(x == min_cor, arr.ind = T)[1,]
  min <- get_names(dn, imin)
  df <- data.frame(imax, imin, max, min)
  return(list(max = max_cor, min = min_cor, df = df))
}


plot_cor_map <- function(x, main = "", panel.title.format = "") {
  get_pts <- function(s, FUN = colMaxs, ...) {
    pts <- rasterToPoints(s, spatial = T)
    nl <- length(names(s))
    df.pts <- as.data.frame(pts)
    m <- FUN(df.pts)
    idx <- sapply(1:nl, function(l) which(df.pts[,l] == m[l], arr.ind = T))
    return(sapply(1:nl, function(l) pts[idx[l], l]))
  }

  xn <- x
  dn  <- dimnames(xn)
  dnm <- names(dn)
  l <- list(); i <- 3
  while (i <= length(dim(xn))) {
    if ( dim(xn)[i] == 1 ) {
      lv <- list(dn[[dnm[i]]][1])
      names(lv) <- dnm[i]
      l <- append(l, lv)
    }
    i <- i + 1
  }
  xn <- drop(xn)
  while (length(dim(xn)) > 3) {
    dn  <- dimnames(xn)
    dnm <- names(dn)
    i <- length(dim(xn))
    lv <- list(dn[[dnm[i]]][1])
    names(lv) <- dnm[i]
    l <- append(l, lv)
    xn <- index_array(xn, i, 1)
  }
  sub <- if (length(l) > 0) paste(names(l), "=", unlist(l), collapse = " | ") else ""
  print(sub)
  if (length(dim(xn)) == 2) xn <- array(xn, dim = c( dim(xn), 1))
  xn <- aperm(xn, c(2,1,3))
  if (!is.character(panel.title.format) || panel.title.format == "")
    panel.title.format = paste(names(dimnames(xn))[3], "= %s")
  panel.title <- sprintf(panel.title.format, dimnames(xn)[[3]])

  lon <- as.numeric(colnames(xn))
  lat <- as.numeric(rownames(xn))
  stck <- stack(apply(xn, 3, raster,
                      xmn = min(lon), xmx = max(lon),
                      ymn = min(lat), ymx =  max(lat)))
  pts_mx <- get_pts(stck, colMaxs)
  pts_mn <- get_pts(stck, colMins)
  levelplot(stck, par.settings = BuRdTheme, names.attr = panel.title, margin = list(NULL),
            xlab = 'Longitude', ylab = 'Latitude', main = main, sub = sub) +
    latticeExtra::layer(sp.polygons(cntry, fill = 'transparent', col = "gray50", alpha = 0.9), under = F) +
    latticeExtra::layer(sp.points(pts_mx[[panel.number()]], pch = 3, cex = 2, lwd = 2, col = "black"), data = list(pts_mx = pts_mx)) +
    latticeExtra::layer(sp.points(pts_mx[[panel.number()]], pch = 1, cex = 2, lwd = 2, col = "black"), data = list(pts_mx = pts_mx)) +
    latticeExtra::layer(sp.points(pts_mn[[panel.number()]], pch = 3, cex = 2, lwd = 2, col = "black"), data = list(pts_mn = pts_mn)) +
    latticeExtra::layer(sp.points(pts_mn[[panel.number()]], pch = 1, cex = 2, lwd = 2, col = "black"), data = list(pts_mn = pts_mn))
}

analyse_cor <- function(x) {
  st <- stat_cor(x)
  cat("Max Correlation = ", st$max, "\n")
  cat("Min Correlation = ", st$min, "\n")
  print(st$df)
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

w   <- load_wind()
pm  <- read_pm10()[, c(4,6)]
hgt <- read_hgt()[,1]
cor_alfa <- get_correlations(w, pm)
wi <- load_interpolated_wind_1000(w)
cor_alfa_wi <- get_correlations(wi, pm)

dom1 <- create_domain(w,
                      lon = c(17.5, 20.0, 22.5,
                              20.0, 22.5, 25.0,
                              22.5, 25.0, 27.5,
                              25.0, 27.5, 30.0),
                      lat = c(32.5, 32.5, 32.5,
                              35.0, 35.0, 35.0,
                              37.5, 37.5, 37.5,
                              40.0, 40.0, 40.0))
m1 <- do_model1(dom1)
# m2 <- do_model2(dom1, 9)
#
dom2 <- create_domain(w, lon = 27.5, lat = 42.5)
m3  <- do_model2(dom2, 9)
# savecsv(x, "1000mb-pm10-log_pm10-hght1-uv-alfa-9.csv")


# plot_vector_field2(u, v, main = "Wind Field")
# plot_streamlines(u, v)


#
# plot_vector(ua, va, dom, m1, level = 1, time = 1, show_grids=T, show_vectors=F)
#
# summary(m1$w$h)
# boxplot(m1$w$h, pch=20, cex=0.3, ylab="speed", xlab="Levels",
#         main="Wind Speed in Domain");abline(h=0, col="red", lty=2)
# hist(m1$w$h[,1], breaks=20, xlab="Vec. Mag.", las=1, col="gray90",
#      main="1000 mb. NVM Histogram")
#
# summary(m1$w$d)
# boxplot(m1$w$d, pch=20, cex=0.3)
# h <- hist(m1$w$d_compass2[,"150"], breaks=seq(0.5,16.5),
#           xlab= "Direction (deg)", las=1, col="lightgreen",
#           main="1000 mb. Wind Direction Histogram", xaxt="n")
# axis(side=1, at=seq(1,16,1), labels=c("N","NNE","NE","ENE",
#                                      "E","ESE", "SE", "SSE",
#                                      "S","SSW","SW","WSW",
#                                      "W","WNW","NW","NNW"))
#
# summary(m2)
# m2[m2 < 0] <- 0
# boxplot(m2, pch=20, cex=0.3, xlab="Levels"); abline(h=0, col="red", lty=2)
# hist_level <- "1000"
# hist(m2[,hist_level], breaks=20, xlab="Southerly", las=1, col="pink",
#      main=paste0(hist_level, " mb. Southerly Wind Histogram"))
#
# # south_wind <- as.vector(which(m1$w$d[,1] >= 124 & m1$w$d[,1] <= 197))
# # south_wind <- as.vector(which(m1$w$d[,1] >=   0 & m1$w$d[,1] <= 361
# #                             & m1$w$h[,1] > 0.95 & m1$w$h[,1] <= 1))
#
# savecsv(m2, "wind_southerly.csv")
# savecsv(m2, "wind_southerly_ma_1.csv")
# savecsv(m1$w$h, "wind_speed.csv")
# savecsv(m1$w$d, "wind_dir.csv")
# savecsv(m1$w$d_compass, "wind_dir_2.csv")
# savecsv(m1$w$d_compass2, "wind_dir_3.csv")
# savecsv(m1$w$d_ma, "wind_dir_moving_average.csv")
# savecsv(m1$w$h_ma, "wind_speed_moving_average.csv")
#
#
# m1wd <- m1$w$d
# m1wdL1 <- m1wd[,1]
#
# m2L1 <- m2[,1]
# south_wind <- which(m2L1 > 0)
#
# summary(m1wdL1[south_wind])
#
# m2SL1 <- m2L1[south_wind]
# m2SL1 <- as.data.frame(m2SL1)
# m2SL1 <- m2SL1[order(m2SL1[,1]),,drop=F]
#
# hist(m2SL1[,1])
# for( i in rownames(tail(m2SL1, 10))) {
#   plot_vector(ua, va, dom, m1, level = 1, time = "2008-01-18")
#   # break;
# }
