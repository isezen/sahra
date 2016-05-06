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
  calc <- function(uv) {
    rot <- function(a) rotax(uv, alfa = a)
    cor2 <- function(x1) {
      ret2 <- apply(pm, 2, cor, y = x1, use = "pairwise.complete.obs")
    }
    r <- vapply(alfa, rot, outer(1:dim(uv)[1], 1:dim(uv)[2]))
    ret <- apply(r, c(2, 3), cor2)
    dimnames(ret)[[3]] <- alfa
    return(ret)
  }
  comp_loc <- which(dim(x) == 2) # find uv dim
  if (length(comp_loc) == 0)
    stop("x does not have a dim = 2 represents u and v")
  time_loc <- which(dim(x) == nrow(pm)) # find time dim
  if (length(time_loc) == 0)
    stop("x does not have a dim suitable with pm")
  if (require(rwind)) {
    if (require(parallel) && par) {
      cl <- parallel::makeCluster(getOption("cl.cores", detectCores()/2))
      parallel::clusterExport(cl, c("pm", "alfa"), envir = environment())
      parallel::clusterEvalQ(cl, library(rwind))
      r <- parallel::parApply(cl, x,
                              (1:length(dim(x)))[-c(comp_loc, time_loc)],
                              calc)
      parallel::stopCluster(cl)
    } else {
      r <- apply(x, (1:length(dim(x)))[-c(comp_loc, time_loc)], calc)
    }
  } else {
    stop("Install rwind package")
  }
  r <- array(r, dim = c(ncol(pm), 2, length(alfa), dim(x)[-c(comp_loc, time_loc)]))
  names(dim(r)) <- c("pm", "component", "alfa", names(dim(x)[-c(comp_loc, time_loc)]))
  dimnames(r) <- append(list(pm = colnames(pm),
                             component = c("u", "v"),
                             alfa = alfa), dimnames(x)[-c(comp_loc, time_loc)])
  dim_order <- 1:length(dim(r))
  comp_loc <- which(dim(r) == 2)
  alfa_loc <- which(dim(r) == length(alfa))
  pm_loc <- which(dim(r) == ncol(pm))
  dim_comp <- dim_order[comp_loc]
  dim_alfa <- dim_order[alfa_loc]
  dim_pm <- dim_order[pm_loc]
  dim_remain <- dim_order[-c(comp_loc, alfa_loc, pm_loc)]
  r <- aperm(r, c(dim_remain, dim_pm, dim_comp, dim_alfa))
  return(r)
}

test_cor2 <- function() {
  nc <- 5
  pm <- matrix(runif(1827 * nc, 0, 150), ncol = nc)
  times <- as.character(seq(from = as.Date("2008-01-01"), length.out = 1827, by = "1 day"))
  rownames(pm) <- times
  colnames(pm) <- paste0("pm", 1:ncol(pm))
  #
  wdim <- c(23, 21, 17, 1827, 2)
  x <- array(runif(prod(wdim), -5, 5), dim = wdim)
  dimnames(x) <- list(lat = 1:wdim[1], lon = 1:wdim[2],
                      level = 1:wdim[3], time = times, comp = c("u", "v"))
  cc <- cor2(x, pm)
}
