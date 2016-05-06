# Saharan Dust Transport Research
# 2016-05-04 Ismail SEZEN
# sezenismail@gmail.com

source("code/base.r")

corstat <- function(x, dim = NULL) {
  csb <- function(x, func = max) {
    dn <- dimnames(x)
    as_char <- function(i) {
      if (is.vector(i)) i <- matrix(i, ncol = 1)
      ret <- sapply(1:ncol(i), function(j) dn[[j]][i[,j]])
      if (is.vector(ret)) ret <- matrix(ret, ncol = length(ret), byrow = T)
      return(ret)
    }
    fn <- as.character(substitute(func))
    fres <- func(x)
    i <-  which(x == fres, arr.ind = T); rownames(i) <- NULL
    ichar <- as_char(i); colnames(ichar) <- paste0('v', colnames(i))
    colnames(i) <- paste0('i', colnames(i))
    df <- data.frame(fres, i, ichar, check.names = F)
    colnames(df)[1] <- fn
    return(df)
  }

  cs <- function(x) {
    mx <- csb(x, max)
    mn <- csb(x, min)
    colnames(mx)[1] <- colnames(mn)[1] <- "max_min"
    df <- rbind(mx, mn)
    return(df)
  }

  if (is.null(dim)) {
    return(cs(x))
  } else {
    l <- list()
    for (i in dimnames(x)[[dim]])
      l[[i]] <- cs(index_array(x, dim, i))
    return(l)
  }
}

corstat_from_files <- function() {
  dir_data_cor <- "data/cor"
  files <- list.files(dir_data_cor, full.names = T)
  for (f in files) {
    bf <- basename(f)
    load(f, environment())
    cat(bf, "\n")
    print(corstat(data))
    cat("\n")
  }
}

cor2 <- function(x, pm, alfa = seq(0, 360, 20), par = T) {
  calc <- function(uv) {
    rot <- function(a) rotax(uv, alfa = a)
    cor2 <- function(x1) {
      apply(pm, 2, cor, y = x1, use = "pairwise.complete.obs")
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
      cl <- parallel::makeCluster(getOption("cl.cores", detectCores() / 2))
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
  r <- array(r, dim = c(ncol(pm), 2,
                        length(alfa),
                        dim(x)[-c(comp_loc, time_loc)]))
  names(dim(r)) <- c("pm",
                     "component",
                     "alfa",
                     names(dim(x)[-c(comp_loc, time_loc)]))
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
  times <- as.character(seq(from = as.Date("2008-01-01"),
    length.out = 1827, by = "1 day"))
  rownames(pm) <- times
  colnames(pm) <- paste0("pm", 1:ncol(pm))
  #
  wdim <- c(23, 21, 17, 1827, 2)
  x <- array(runif(prod(wdim), -5, 5), dim = wdim)
  dimnames(x) <- list(lat = 1:wdim[1], lon = 1:wdim[2],
                      level = 1:wdim[3], time = times, comp = c("u", "v"))
  cc <- cor2(x, pm)
}
