# Saharan Dust Transport Research
# 2016-05-04 Ismail SEZEN
# sezenismail@gmail.com

source("code/correlation.r")
source("code/filehelper.r")


calcor <- function(file = stop("'file' must be specified")) {
  w <- load_rdata(file)
  pm <- read_pm10()
  for (i in 1:ncol(pm)) {
    p <- pm[,i]
    pname <- colnames(p)
    correlation <- cor2(w, p, alfa = seq(0, 360, 1))
    save_rdata(correlation, file = paste0("cor_", pname, "_", basename(file)))
  }
}

calcor(file = "data/rdata/r2-pres-00-uv.rdata")
calcor(file = "data/rdata/r2-pres-06-uv.rdata")
calcor(file = "data/rdata/r2-pres-12-uv.rdata")
calcor(file = "data/rdata/r2-pres-18-uv.rdata")
calcor(file = "data/rdata/r2-pres-daily-uv.rdata")
calcor(file = "data/rdata/r2-pres-daily-hgt.rdata")
