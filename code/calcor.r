# Saharan Dust Transport Research
# 2016-05-04 Ismail SEZEN
# sezenismail@gmail.com

source("code/correlation.r")
source("code/filehelper.r")

calcor <- function(files = stop("'file' must be specified")) {
  pm <- read_pm10()
  dir_out <- "data/cor"
  dir.create(dir_out, showWarnings = F)
  nof <- length(files)
  i <- 1
  for (f in files) {
    w <- load_rdata(f)
    fwe <- basename(tools::file_path_sans_ext(f))
    save_to <- file.path(dir_out, paste0("cor_", fwe, ".rds"))
    if (!file.exists(save_to)) {
      cat("Calculation started at", as.character(now()), "\n")
      cat("(", i, "/", nof, ") ", basename(save_to), "\n", sep = "")
      st <- system.time(data <- cor2(w, pm, alfa = seq(0, 360, 1)))[3]
      cat("[Elapsed :", st, "sec]\n")
      saveRDS(data, file = save_to)
    } else {
      cat("(", i, "/", nof, ") ", save_to, " is exist.\n", sep = "")
    }
    i <- i + 1
  }
}

files <- c("r2-pres-00-uv.rdata", "r2-pres-06-uv.rdata",
           "r2-pres-12-uv.rdata", "r2-pres-18-uv.rdata",
           "r2-pres-daily-uv.rdata",
           "r2-10m-00-uv.rdata", "r2-10m-06-uv.rdata",
           "r2-10m-12-uv.rdata", "r2-10m-18-uv.rdata",
           "r2-10m-daily-uv.rdata",
           "r2-pres-00-hgt.rdata", "r2-pres-06-hgt.rdata",
           "r2-pres-12-hgt.rdata", "r2-pres-18-hgt.rdata",
           "r2-pres-daily-hgt.rdata")
files <- file.path("data/rds", files)
calcor(files)
