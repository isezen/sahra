# Saharan Dust Transport Research
# 2016-05-04 Ismail SEZEN
# sezenismail@gmail.com

source("code/correlation.r")

calcor <- function(files = stop("'file' must be specified")) {
  pm <- read_pm10()
  dir_out <- "data/cor"
  dir.create(dir_out, showWarnings = F)
  nof <- length(files)
  i <- 1
  for (f in files) {
    fwe <- basename(tools::file_path_sans_ext(f))
    save_to <- file.path(dir_out, paste0("cor_", fwe, ".rds"))
    if (!file.exists(save_to)) {
      cat("Calculation started at", as.character(now()), "\n")
      cat("(", i, "/", nof, ") ", basename(save_to), "\n", sep = "")
      w <- readRDS(f)
      st <- system.time(data <- cor2(w, pm, alfa = seq(0, 360, 1)))[3]
      cat("[Elapsed :", st, "sec]\n")
      saveRDS(data, file = save_to)
    } else {
      cat("(", i, "/", nof, ") ", save_to, " is exist.\n", sep = "")
    }
    i <- i + 1
  }
}

files <- paste0(c("r2-pres-00-uv", "r2-pres-06-uv",
                  "r2-pres-12-uv", "r2-pres-18-uv",
                  "r2-pres-daily-uv",
                  "r2-surf-00-uv", "r2-surf-06-uv",
                  "r2-surf-12-uv", "r2-surf-18-uv",
                  "r2-surf-daily-uv",
                  "r2-pres-00-hgt", "r2-pres-06-hgt",
                  "r2-pres-12-hgt", "r2-pres-18-hgt",
                  "r2-pres-daily-hgt",
                  "r2-surf-00-mslp", "r2-surf-06-mslp",
                  "r2-surf-12-mslp", "r2-pres-18-mslp",
                  "r2-surf-daily-mslp",
                  "r2-pres-00-omega", "r2-pres-06-omega",
                  "r2-pres-12-omega", "r2-pres-18-omega",
                  "r2-pres-daily-omega"), ".rds")
files <- file.path("data/rds", files)
calcor(files)
