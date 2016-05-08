# Saharan Dust Transport Research
# 2016-05-04 Ismail SEZEN
# sezenismail@gmail.com

library(ncdf4)
library(ncdf4.helpers)
library(PCICt)

nc_get_data <- function(nc, lat = NULL, lon = NULL, level = NULL, time = NULL) {
  c_to_i <- function(dn, vals) match(vals, ncvar_get(nc, dn))
  i_to_c <- function(dn, indices) ncvar_get(nc, dn)[indices]
  var_name <- nc$var[[length(nc$var)]]$name

  i <- l <- list()
  i$lon <- if (is.null(lon)) 1:nc$dim$lon$len else c_to_i("lon", lon)
  l$lon  <- i_to_c("lon", i$lon)

  i$lat <- if (is.null(lat)) 1:nc$dim$lat$len else c_to_i("lat", lat)
  l$lat  <- i_to_c("lat", i$lat)

  if ("level" %in% names(nc$dim)) {
    if (is.null(level)) i$level <- 1:nc$dim$level$len
    if (is.character(level)) i$level <- c_to_i("level", level)
    l$level <- i_to_c("level", i$level)
  }

  ts <- nc.get.time.series(nc)
  if (is.character(time)) i$time <- which(as.character(ts) == time)
  if (is.null(time)) i$time <- 1:nc$dim$time$len
  l$time <- as.character(ts[i$time])

  data <- nc.get.var.subset.by.axes(nc, var_name, i,
                                    nc.get.dim.names(nc)[1:length(l)])
  dimnames(data) <- l
  attr(data, 'name') <- var_name
  attr(data, 'longname') <- nc$var[[length(nc$var)]]$longname
  attr(data, 'unit') <- nc$var[[length(nc$var)]]$unit
  attr(data, 'filename') <- nc$filename
  return(data)
}

nc2rds <- function(file_prefix = "r2-pres-4-", suffix = c("u", "v")) {
  l <- list()
  dir_out <- "data/rds"
  dir.create(dir_out, showWarnings = F)
  save_to <- file.path(dir_out, paste0(file_prefix,
                                       paste0(suffix, collapse = ""), ".rds"))
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
      attr(data, "name") <- paste0(suffix, collapse = "")
      attr(data, 'longname') <- gsub("U-wind", "Wind",
                                     attr(l[[1]], 'longname'))
      attr(data, 'unit') <- attr(l[[1]], 'unit')
      attr(data, 'filename') <- sapply(l, attr, "filename")
    } else {
      data <- l[[1]]
    }
    saveRDS(data, save_to)
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
  attrs <- names(attributes(x))
  attrs <- attrs[!(attrs %in% c("dim", "dimnames"))]
  for (a in attrs) attr(e, a) <- attr(x, a)
  attr(e, 'longname') <- gsub("6-hourly ", "", attr(e, 'longname'))
  tm <- substr(dn$time[1], nchar(dn$time[1]) - 7, nchar(dn$time[1]) - 6)
  attr(e, 'longname') <- paste0(attr(e, 'longname'), " at ", tm, " UTC")
  return(e)
}

nc2rds_all <- function() {
  for (f in c("r2-surf-4-", "r2-surf-daily-",
              "r2-pres-4-", "r2-pres-daily-")) nc2rds(f)

  for (f in c("r2-pres-4-", "r2-pres-daily-"))
    for (s in c("hgt", "omega")) nc2rds(f, s)

  for (f in c("r2-surf-4-", "r2-surf-daily-")) nc2rds(f, "mslp")

  dir_out <- "data/rds"
  for (f in list.files(dir_out, "*-4-*")) {
    x <- readRDS(file.path(dir_out, f))
    tm <- c("00", "06", "12", "18")
    for (i in 1:4) {
      data <- extract_time(x, from = i)
      save_to <- file.path(dir_out, gsub("-4-", paste0("-", tm[i], "-"), f))
      attr(data, 'filename') <- save_to
      saveRDS(data, save_to)
    }
  }
}
