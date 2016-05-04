# Saharan Dust Transport Research
# 2016-05-04 Ismail SEZEN
# sezenismail@gmail.com

windhelper_get_args <- function(...) {
  l <- list(...); x <- NULL
  if (length(l) == 1) {
    if (is.array(l[[1]])) x <- l[[1]]
    if (is.data.frame(l[[1]])) x <- as.matrix(l[[1]])
  }
  if (is.null(x)) x <- do.call(c, l)
  if (is.vector(x) && (length(x) %% 2) == 0)
    x <- array(x, dim = c(length(x) / 2, 2))
  return(x)
}


deg_to_comp <- function(x, bin_names=c("N", "NNE", "NE", "ENE",
                                       "E", "ESE", "SE", "SSE",
                                       "S", "SSW", "SW", "WSW",
                                       "W", "WNW", "NW", "NNW")) {
  x[x == 0] <- 360
  i <- 360 / length(bin_names)
  val <- round( (x / i) + .5)
  return(matrix(bin_names[val], nrow(x), ncol(x), dimnames = dimnames(x)))
}

# Convert u,v wind components to horizontal speed and direction
uv2hd <- function(..., degree = T, reverse = T, drop = T) {
  source("code/base.r", local = T)
  x <- windhelper_get_args(...)
  comp_loc <- which(dim(x) == 2, arr.ind = T) # find uv dim order
  if (length(comp_loc) == 0)
    stop("x does not have a dim = 2 represents u and v")
  dx <- dim(x); ldx <- length(dx)
  r <- aperm(x, c(comp_loc, (1:ldx)[-comp_loc])) # get uv dim to start

  num_pi <- if (degree) 180 else pi
  add_pi <- if (reverse) num_pi else 0
  u <- index_array(r, 1, 1)
  v <- index_array(r, 1, 2)
  h <- sqrt(colSums(r ^ 2))
  d <- ((atan2(u / h, v / h) * num_pi / pi) + add_pi)
  if (degree) d <- d  %% 360
  r <- array(c(h, d), dim = dx)

  r[is.nan(r)] <- 0
  dmn <- dimnames(x)
  if (is.null(dmn)) {
    dimnames(r) <- list(NULL, c("h", "d"))
  } else {
    dimnames(r) <- dmn
  }
  if (drop) r <- drop(r)
  return(r)
}

hd2uv <- function(..., degree = T, drop = T) {
  source("code/base.r", local = T)
  x <- windhelper_get_args(...)
  comp_loc <- which(dim(x) == 2, arr.ind = T) # find uv dim order
  if (length(comp_loc) == 0)
    stop("x does not have a dim = 2 represents h and d")
  dx <- dim(x); ldx <- length(dx)
  r <- aperm(x, c(comp_loc, (1:ldx)[-comp_loc])) # get uv dim to start

  num_pi <- if (degree) 180 else pi
  h <- index_array(r, 1, 1)
  d <- index_array(r, 1, 2)
  d <- (d - num_pi) / num_pi
  u <- h * sinpi(d)
  v <- h * cospi(d)
  r <- array(c(u, v), dim = dx)

  r[is.nan(r)] <- 0
  dmn <- dimnames(x)
  if (is.null(dmn)) {
    dimnames(r) <- list(NULL, c("h", "d"))
  } else {
    dimnames(r) <- dmn
  }
  if (drop) r <- drop(r)
  return(r)
}

rotax <- function(..., alfa=45, degree=T, right=T, drop = T) {
  source("code/base.r", local = T)
  x <- windhelper_get_args(...)
  comp_loc <- which(dim(x) == 2, arr.ind = T) # find uv dim order
  if (length(comp_loc) == 0)
    stop("x does not have a dim = 2 represents u and v")

  u <- index_array(x, comp_loc, 1)
  v <- index_array(x, comp_loc, 2)
  a <- if (degree) alfa / 180 else alfa / pi
  c <- cospi(a); s <- sinpi(a)
  if (right) {
    un <- u * c - v * s
    vn <- u * s + v * c
  } else {
    un <- u * c + v * s
    vn <- v * c - u * s
  }
  dm <- dim(x); dn <- dimnames(x)
  x <- array(c(un, vn), dim = dm)
  dimnames(x) <- dn
  if (drop) x <- drop(x)
  return(x)
}
