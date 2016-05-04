# Saharan Dust Transport Research
# 2016-05-04 Ismail SEZEN
# sezenismail@gmail.com


ibicubic <- function(x, xlim = NULL, ylim = xlim, dx = NULL, dy = dx, par = T) {
  stopifnot( length(dim(x)) > 1 )

  Rev <- function(x, margin, ...) {
    if (!is.array(x))
      stop("'x' is not an array")
    newdim <- rep("", length(dim(x)))
    newdim[margin] <- paste(dim(x), ":1", sep = "")[margin]
    z <- eval(parse(text = gettextf("x[%s]", paste(newdim, sep = "",
                                                   collapse = ","))))
    class(z) <- oldClass(x)
    return(z)
  }

  set_dim_monoton_inc <- function(x) {
    dims <- vector()
    for (i in 1:2) { # check only first two dims
      n <- as.numeric(dimnames(x)[[i]])
      if (length(n)) {
        x <- Rev(x, i)
        dims <- c(dims, i)
      } else {
        dimnames(x)[i] <- list(1:dim(x)[i])
      }
    }
    return(list(x = x, dims = dims))
  }

  new_xy <- function(x, xlim, dx) {
    rng <- range(x)
    new_r <- seq(rng[1], rng[2], dx)
    return(new_r[new_r >= xlim[1] & new_r <= xlim[2]])
  }

  x <- set_dim_monoton_inc(x)
  dims <- x$dims; x <- x$x

  r <- as.numeric(dimnames(x)[[1]])
  c <- as.numeric(dimnames(x)[[2]])
  if (is.null(dx)) dx <- abs((r[2] - r[1]) / 2)
  if (is.null(dy)) dy <- abs((c[2] - c[1]) / 2)
  if (is.null(xlim)) xlim <- range(r)
  if (is.null(ylim)) ylim <- range(c)
  new_r <- new_xy(r, xlim, dx)
  new_c <- new_xy(c, ylim, dy)
  l <- length(dim(x))
  if (l > 2) {
    if (par) {
      library(parallel)
      cl <- makeCluster(detectCores() - 1)
      clusterExport(cl, c("r", "c", "xlim", "ylim", "dx", "dy"), envir = environment())
      clusterEvalQ(cl, library(akima))
      xi <- parApply(cl, x, (1:l)[-c(1,2)], function(z) bicubic.grid(r, c, z, xlim, ylim, dx, dy)$z)
      stopCluster(cl)
    } else {
      library(akima)
      xi <- apply(x, (1:l)[-c(1,2)], function(z) bicubic.grid(r, c, z, xlim, ylim, dx, dy)$z)
      detach('package:akima')
    }
  } else {
    library(akima)
    xi <- bicubic.grid(r, c, x, xlim, ylim, dx, dy)$z
    detach('package:akima')
  }
  xi <- array(xi, dim = c(length(new_r), length(new_c), dim(x)[-c(1,2)]))
  dimnames(xi)[[1]] <- new_r
  dimnames(xi)[[2]] <- new_c
  dimnames(xi)[-c(1,2)] <- dimnames(x)[-c(1,2)]
  if (length(dims)) xi <- Rev(xi, dims)
  return(xi)
}

