# Saharan Dust Transport Research
# 2016-05-04 Ismail SEZEN
# sezenismail@gmail.com

library(ncdf4)
library(ncdf4.helpers)

source("code/filehelper.r")

nc_get_data <- function(nc, lats=NULL, lons=NULL, levels=NULL, times=NULL) {
  c_to_i <- function(dim.name, vals) {
    return(match(vals, ncvar_get(nc, dim.name)))
  }
  i_to_c <- function(dimname, indices) {
    return(ncvar_get(nc, dimname)[indices])
  }

  var_name <- nc$var[[length(nc$var)]]$name

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

nc2rdata <- function(file_prefix = "r2-pres-4-", suffix = c("u", "v")) {
  l <- list()
  for (c in suffix) {
    f <- paste0(file_prefix, c, ".nc")
    nc <- nc_open(file.path("data", "nc", f))
    l[[c]]  <- nc_get_data(nc)
    nc_close(nc)
  }
  if (length(l) > 1) {
    data <- array(unlist(l), dim = c(dim(l[[1]]), component = length(l)))
    dimnames(data) <- append(dimnames(l[[1]]), list(component = suffix))
  } else {
    data <- l[[1]]
  }
  save_rdata(data, file = paste0("rdata/", file_prefix, paste0(suffix, collapse = "")))
}
