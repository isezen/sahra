# Saharan Dust Transport Research
# 2016-05-04 Ismail SEZEN
# sezenismail@gmail.com

download_correlations <- function() {
  dir_out <- "data/cor"
  dir.create(dir_out, showWarnings = F)
  system2('scp', paste0('mustang:/home/isezen/proj/sahra/data/cor/\\*.rds ',
                        dir_out))
}

download_ncep_R2 <- function() {
  url <- "https://kovan.itu.edu.tr/index.php/s/qCe6lasJdaf3MxC/download"
  file <- "sahra-ncep-r2.zip"
  download.file(url, file, method = "auto")
  fl <- unzip(file, list = T)
  fl <- fl$Name
  unzip(file, fl[!grepl("__MACOSX/", fl)], exdir = "data/nc")
}

save_rdata <- function(..., file = stop("'file' must be specified")) {
  if (nchar(tools::file_ext(file)) == 0) file <- paste0(file, ".rdata")
  save(..., file = file.path("data/rdata", file), envir = parent.frame())
}

load_rdata <- function(file = stop("'file' must be specified")) {
  load(file)
  get(ls()[ls() != file])
}

save_csv <- function(x, file = stop("'file' must be specified"),
                     format = "%d.%m.%Y") {
  library(xts)
  dts <- if (is.xts(x)) index(x) else as.POSIXct(rownames(x), tz = "GMT")
  dts <- format(dts, format)
  df <- as.data.frame(x)
  df <- cbind(dts, df)
  colnames(df)[[1]] <- "date"
  write.table(df, file = file, quote = F, sep = ";", row.names = F)
}

read_csv <- function(file = stop("'file' must be specified"), shift = 0, filter = NULL) {
  d <- read.table(file.path("data/csv", file))
  # find PM
  if (!is.null(filter)) {
    indices <- c(5, which(grepl(filter, colnames(d)) == T))
    d <- d[, indices] # filter only PM data
    d <- d[, colSums(is.na(d)) < nrow(d)] # remove full NA cols
  }
  library(lubridate)
  library(zoo)
  library(xts)
  d[, 1] <- as.POSIXct(d[, 1], tz = "GMT") + days(shift)
  d <- xts(d[, -1], order.by = d[, 1])

  dts <- seq(as.POSIXct("2008-01-01", tz = "GMT"),
             as.POSIXct("2015-12-27", tz = "GMT"), by = "days")
  d <- merge(d, zoo(, dts), all = T)
  d <- d["2008/2012"]
  return(d)
}

read_pm10 <- function() read_csv("NightResults12042016.csv", 0, "PM")
read_wind <- function() read_csv("NightResults12042016.csv", 0, "Wind")
read_hgt <- function() read_csv("NightResults12042016.csv", 0, "HGHT")
