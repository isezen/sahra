# Saharan Dust Transport Research
# 2016-05-04 Ismail SEZEN
# sezenismail@gmail.com

source("code/correlation.r")
source("code/filehelper.r")


calcor <- function(files = stop("'file' must be specified")) {
  pm <- read_pm10()
  dir_out <- "data/cor"
  dir.create(dir_out, showWarnings = F)
  for (f in files) {
    w <- load_rdata(f)
    save_to <- file.path(dir_out, paste0("cor_", basename(f)))
    if (!file.exists(save_to)) {
      cat("Calculation started at ", as.character(now()), "\n")
      cat(basename(save_to), "\n")
      print(system.time(data <- cor2(w, pm, alfa = seq(0, 360, 20))))
      save(data, file = save_to, envir = environment())
    } else {
      cat(save_to, " is exist.\n")
    }
  }
}

files <- c("r2-pres-00-uv.rdata", "r2-pres-06-uv.rdata",
           "r2-pres-12-uv.rdata", "r2-pres-18-uv.rdata",
           "r2-pres-daily-uv.rdata",
           "r2-10m-00-uv.rdata", "r2-10m-06-uv.rdata",
           "r2-10m-12-uv.rdata", "r2-10m-18-uv.rdata",
           "r2-10m-daily-uv.rdata")
files <- file.path("data/rdata", files)
calcor(files)
