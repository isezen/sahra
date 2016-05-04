# Saharan Dust Transport Research
# 2016-05-04 Ismail SEZEN
# sezenismail@gmail.com

load_countries <- function() {
  # download map
  dir_countries <- "countries_shp"
  fl <- "ne_50m_admin_0_countries"
  if (!dir.exists(dir_countries)) {
    fle <- paste0(fl, ".zip")
    url <- paste0("http://naciscdn.org/naturalearth/50m/cultural/", fle)
    download.file(url, fle, method = "curl")
    unzip(fle, exdir = dir_countries)
  }
  library(rgdal)
  cntry <- readOGR(dsn = paste0("./", dir_countries), verbose = F,
                   layer = "ne_50m_admin_0_countries", stringsAsFactors = T)
  return(cntry)
}
cntry <- load_countries()


plot_cor_map <- function(x, main = "", panel.title.format = "") {
  library(raster)
  library(rasterVis)
  library(timeSeries)
  source("code/base.r", local = T)
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

plot_vector_field <- function(u, v, ...) {
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
