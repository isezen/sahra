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
  library(latticeExtra)
  library(raster)
  library(rasterVis)
  library(timeSeries)
  get_pts <- function(s, FUN = colMaxs, ...) {
    pts <- rasterToPoints(s, spatial = T)
    nl <- length(names(s))
    df.pts <- as.data.frame(pts)
    m <- FUN(df.pts)
    idx <- sapply(1:nl, function(l) which(df.pts[,l] == m[l], arr.ind = T))
    return(sapply(1:nl, function(l) pts[idx[l], l]))
  }

  xn <- rsezen::ibicubic(x, dx = 0.05)
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
    xn <- rwind:::index_array(xn, i, 1)
  }
  sub <- if (length(l) > 0) paste(names(l), "=", unlist(l), collapse = " | ") else ""
  # print(sub)
  if (length(dim(xn)) == 2) {
    xn <- array(xn, dim = c( dim(xn), 1))
    names(dim(xn)) <- dnm[1:3]
    dimnames(xn) <- dn[1:3]
  }
  xn <- aperm(xn, c(2,1,3))
  if (!is.character(panel.title.format) || panel.title.format == "")
    panel.title.format = paste(names(dimnames(xn))[3], "= %s")
  panel.title <- sprintf(panel.title.format, dimnames(xn)[[3]])

  lon <- as.numeric(colnames(xn))
  lat <- as.numeric(rownames(xn))
  # print(lon); print(lat)
  stck <- stack(apply(xn, 3, raster,
                      xmn = min(lon), xmx = max(lon),
                      ymn = min(lat), ymx =  max(lat)))
  pts_mx <- get_pts(stck, colMaxs)
  pts_mn <- get_pts(stck, colMins)
  levelplot(stck, par.settings = BuRdTheme, names.attr = panel.title, margin = list(NULL),
            xlab = 'Longitude', ylab = 'Latitude', main = main, sub = list(label = sub, cex = 0.7)) +
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

plot_cor_for <- function(pattern = "hgt") {
  dir_data_cor <- "data/cor"
  files <- list.files(dir_data_cor, full.names = T)
  for (p in pattern) files <- Filter(function(x) grepl(p, x), files)
  for (f in files) {
    bf <- basename(f)
    if (tools::file_ext(bf) == "rds") {
      x <- readRDS(f)
      df <- corstat(x)
      cap <- attr(x, "longname")
      if (is.null(cap)) cap <- bf
      panel.title.format <- ""
      imx <- which(df[,1] == max(df[,1]), arr.ind = T)[1]
      df <- df[imx,]
      if (grepl("uv", bf)) {
        if (grepl("pres", bf)) {
          x <- x[,, df[1, "ipm"],, df[1, "ialfa"], df[1, "ilevel"], drop = F]
        } else if (grepl("surf", bf)) {
          x <- x[,, df[1, "ipm"],, df[1, "ialfa"], drop = F]
        }
        panel.title.format <- "%s Wind"
      } else if (grepl("hgt", bf)) {
        x <- x[,,,  df[1, "ilevel"], drop = T]
        x <- x[,,df[1, "ipm"], drop = F]
      } else if (grepl("mslp", bf)) {
        x <- x[,, df[1, "ipm"], drop = F]
      } else if (grepl("omega", bf)) {
        x <- x[,,,  df[1, "ilevel"], drop = T]
        x <- x[,,df[1, "ipm"], drop = F]
      }
      pdf(paste0(bf, ".pdf"), width = 16, height = 8)
      print(plot_cor_map(x, main = cap, panel.title.format = panel.title.format))
      dev.off()
    }
  }
}

plot_reg <- function(pattern = NULL, save = T) {
  dir_data_cor <- "data/cor"
  files <- list.files(dir_data_cor, full.names = T)
  for (p in pattern) files <- Filter(function(x) grepl(p, x), files)
  pm <- read_pm10()
  for (f in files) {
    bf <- basename(f)
    if (tools::file_ext(bf) == "rds") {
      tryCatch({
        xc <- readRDS(f)
        longname <- attr(xc, "longname")
        varname <- attr(xc, "name")
        if (is.null(longname)) longname <- bf
        df <- corstat(xc)
        imx <- which(df[,1] == max(df[,1]), arr.ind = T)[1]
        df <- df[imx,]
        x <- readRDS(attr(xc, "filename"))
        sub <- paste0("lon = ", df$vlon[1], ", lat = ", df$vlat[1])
        if (grepl("uv", bf)) {
          sub <- paste0(sub, ", alfa = ", df$valfa[1])
          x <- rwind::rotax(x, alfa = as.numeric(as.character(df$valfa[1])))
          if (grepl("pres", bf)) {
            xlab <- paste0(df$vcomponent[1], "-", longname)
            xlab <- gsub("Pressure Levels", paste0(df$vlevel, " hPa"), xlab)
            x <- x[df$ilon[1], df$ilat[1], df$ilevel[1],,df$icomponent[1]]
            xc <- xc[,, df$ipm[1], df$icomponent[1], df$ialfa[1], df$ilevel[1], drop = F]
          } else if (grepl("surf", bf)) {
            if (grepl("6-Hourly Forecast of ", longname)) {
              xlab <- gsub("6-Hourly Forecast of ", paste0(df$vcomponent[1], "-"), longname)
            } else {
              xlab <- gsub("Daily Forecast of ", paste0(df$vcomponent[1], "-"), longname)
            }
            x <- x[df$ilon[1], df$ilat[1],,df$icomponent[1]]
            xc <- xc[,, df$ipm[1], df$icomponent[1], df$ialfa[1], drop = F]
          }
        } else if (grepl("hgt", bf)) {
          xlab <- gsub("Pressure Levels", paste0(df$vlevel, " hPa"), longname)
          x <- x[df$ilon[1], df$ilat[1], df$ilevel[1],]
          xc <- xc[,,, df$ilevel[1], drop = T]
          xc <- xc[,,df$ipm[1], drop = F]
        } else if (grepl("mslp", bf)) {
          x <- x[df$ilon[1], df$ilat[1],] / 100
          xc <- xc[,, df$ipm[1], drop = F]
          xlab <- paste0(longname, " (hPa)")
        } else if (grepl("omega", bf)) {
          xlab <- paste0(longname, " (", attr(x, "unit"), ")")
          x <- x[df$ilon[1], df$ilat[1], df$ilevel[1],]
          xc <- xc[,,, df$ilevel[1], drop = T]
          xc <- xc[,,df$ipm[1], drop = F]
        }
        p <- as.vector(pm[,df$ipm[1]])
        ylab <- as.character(df$vpm[1])
        clr <- "cornflowerblue"
        if  (grepl("log", bf)) {
          p <- log(p)
          ylab <- paste0("log(", ylab, ")")
          clr <- "pink"
        }
        reg <- lm(p ~ x)
        cr <- cor(p, x, use = "pairwise.complete.obs")
        sr <- summary(reg)

        if (save) pdf(paste0("reg_", bf, ".pdf"), width = 12, height = 15)
        def.par <- par(no.readonly = TRUE)
        layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow = TRUE))
        plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
        plot(reg$model[,2:1], pch = 21, bg = clr, las = 1,
             ylab = ylab, xlab = xlab, sub = sub)
        abline(reg, col = "red")
        cf <- round(coef(reg), 5)
        comp_name <- as.character(if (is.null(df$vcomponent[1])) varname else df$vcomponent[1])
        eq <- paste0("pm = ", cf[1],
                     ifelse(sign(cf[2]) == 1, " + ", " - "), abs(cf[2]), " ", comp_name)
        adj.rs <- round(sr$adj.r.squared, 4)
        cr <- round(cr, 3)
        legend("topleft",
               c(eq, paste0("Adj. R^2 = ", adj.rs), paste0("cor = ", cr)),
               bty = "n", col = 1:3, pch = 16, text.col = 1:3, y.intersp = 2)
        plot(reg, ask = F, pch = 21, bg = "chartreuse3")
        print(plot_cor_map(xc, main = "", panel.title.format = ""), split = c(1, 1, 2, 3), newpage = F)
        par(def.par)
        if (save) dev.off()
      }, error = function(e) {
        cat(bf, "cannot be read.\n")
      }, finally = {
        next
      })
    }
  }
}
