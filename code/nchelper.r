# Saharan Dust Transport Research
# 2016-05-04 Ismail SEZEN
# sezenismail@gmail.com

library(ncdf4)
library(ncdf4.helpers)

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
  dir_out <- "data/rdata"
  dir.create(dir_out, showWarnings = F)
  save_to <- file.path(dir_out, paste0(file_prefix, paste0(suffix, collapse = ""), ".rdata"))
  if (!file.exists(save_to)) {
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
    save_to <- file.path(dir_out, paste0(file_prefix, paste0(suffix, collapse = ""), ".rdata"))
    save(data, file = save_to, envir = environment())
  } else {
    cat(save_to, "is exist.\n")
  }
}

extract_time <- function(x, from = 1, by = 4) {
  source("code/base.r", local = T)
  dim_time <- which(names(dim(x)) == "time")
  time_length <- dim(x)[dim_time]
  e <- index_array(x, dim_time, seq(from, time_length, by))
  dn <- dimnames(e)
  names(dim(e)) <- names(dn)
  dimnames(e) <- dn
  return(e)
}

nc2rdata_all <- function() {
  dir_out <- "data/rdata"
  for (f in c("r2-10m-4-", "r2-10m-daily-",
              "r2-pres-4-", "r2-pres-daily-"))
    nc2rdata(f)
  for (f in c("r2-pres-4-", "r2-pres-daily-")) nc2rdata(f, "hgt")
  for (f in list.files(dir_out, "*-4-*")) {
    load(file.path(dir_out, f))
    x <- data
    tm <- c("00", "06", "12", "18")
    for (i in 1:4) {
      data <- extract_time(x, from = i)
      save_to <- file.path(dir_out, gsub("-4-", paste0("-", tm[i], "-"), f))
      save(data, file = save_to, envir = environment())
    }
  }
}
