# Saharan Dust Transport Research
# 2016-05-04 Ismail SEZEN
# sezenismail@gmail.com

corstat <- function(x) {
  get_names <- function(dn, i) {
    m <- vector()
    for (j in 1:length(i)) m <- c(m, dn[[names(i[j])]][i[j]])
    return(m)
  }
  dn <- dimnames(x)
  max_cor <- max(x); min_cor <- min(x)
  imax <-  which(x == max_cor, arr.ind = T)[1,]
  max <- get_names(dn, imax)
  imin <- which(x == min_cor, arr.ind = T)[1,]
  min <- get_names(dn, imin)
  df <- data.frame(imax, imin, max, min)
  return(list(max = max_cor, min = min_cor, df = df))
}

#' Shows a summary of correlation data.
#'
#' @param x correlation array
#' @return Nothing
#' @examples
#' x <- array(rnorm(1000), dim = c(10, 10, 10))
#' cor_anal(x)
cor_anal <- function(x) {
  st <- corstat(x)
  cat("Max Correlation = ", st$max, "\n")
  cat("Min Correlation = ", st$min, "\n")
  print(st$df)
}

cor_anal_from_files <- function() {
  dir_data_cor <- "data/cor"
  files <- list.files(dir_data_cor, full.names = T)
  for (f in files) {
    bf <- basename(f)
    load(f, environment())
    cat(bf, "\n")
    cor_anal(data)
    cat("\n")
  }
}

cor2 <- function(x, pm, alfa = seq(0, 360, 20), par = T) {
  source("code/base.r", local = T)
  source("code/windhelper.r", local = T)
  calc <- function(uv) {
    rot <- function(a) rotax(uv, alfa = a)
    cor2 <- function(x, y) cor(y, x, use = "pairwise.complete.obs")
    r <- vapply(alfa, rot, outer(1:dim(uv)[1], 1:dim(uv)[2]))
    ret <- t(apply(r, c(2, 3), cor, y = pm, use = "pairwise.complete.obs"))
    return(ret)
  }
  comp_loc <- which(dim(x) == 2) # find uv dim order
  if (length(comp_loc) == 0)
    stop("x does not have a dim = 2 represents u and v")
  time_loc <- which(dim(x) == nrow(pm)) # find time dim order
  if (length(time_loc) == 0)
    stop("x does not have a dim suitable with pm")
  if (require(parallel) && par) {
    cl <- makeCluster(getOption("cl.cores", detectCores()/2))
    clusterExport(cl, c("pm"), envir = environment())
    clusterEvalQ(cl, library(rwind))
    r <- parApply(cl, x, (1:length(dim(x)))[-c(comp_loc, time_loc)], calc)
    stopCluster(cl)
  } else {
    r <- apply(x, (1:length(dim(x)))[-c(comp_loc, time_loc)], calc)
  }
  r <- array(r, dim = c(length(alfa), 2, dim(x)[-c(comp_loc, time_loc)]))
  names(dim(r)) <- c("alfa", "component", names(dim(x)[-c(comp_loc, time_loc)]))
  dimnames(r) <- append(list(alfa = alfa, component = c("u", "v")), dimnames(x)[-c(comp_loc, time_loc)])
  dim_order <- 1:length(dim(r))
  comp_loc <- which(dim(r) == 2)
  alfa_loc <- which(dim(r) == length(alfa))
  dim_comp <- dim_order[comp_loc]
  dim_alfa <- dim_order[alfa_loc]
  dim_remain <- dim_order[-c(comp_loc, alfa_loc)]
  r <- aperm(r, c(dim_remain, dim_comp, dim_alfa))
  return(r)
}
