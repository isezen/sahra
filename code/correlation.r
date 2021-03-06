# Saharan Dust Transport Research
# 2016-05-04 Ismail SEZEN
# sezenismail@gmail.com

source("code/filehelper.r")

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
      l[[i]] <- cs(rwind:::index_array(x, dim, i))
    return(l)
  }
}

corstat_from_files <- function(pattern = NULL) {
  dir_data_cor <- "data/cor"
  files <- list.files(dir_data_cor, full.names = T)
  for (p in pattern) files <- Filter(function(x) grepl(p, x), files)
  for (f in files) {
    bf <- basename(f)
    if (tools::file_ext(bf) == "rds") {
      tryCatch({
        data <- readRDS(f)
        longname <- attr(data, "longname")
        if (is.null(longname)) longname <- bf
        if (grepl("log", bf)) longname <- paste0(longname, " (log)")
        cat(longname, ":\n")
        print(corstat(data))
        cat("\n")
      }, error = function(e) {
        cat(bf, "cannot be read.\n")
      }, finally = {
        next
      })
    }
  }
}

get_data_by_cor <- function(fun=max) {
  dir_data_cor <- "data/cor"
  files <- list.files(dir_data_cor, full.names = T)
  files <- Filter(function(x) grepl("daily", x), files)
  ret <- NULL
  for (f in files) {
    bf <- basename(f)
    if (tools::file_ext(bf) == "rds") {
      tryCatch({
        xc <- readRDS(f)
        longname <- attr(xc, "longname")
        varname <- attr(xc, "name")
        if (is.null(longname)) longname <- bf
        if (grepl("log", bf)) longname <- paste0(longname, " (log)")
        cat(longname, ":\n")
        df <- corstat(xc)
        imx <- which(df[,1] == fun(df[,1]), arr.ind = T)[1]
        df <- df[imx,]
        x <- readRDS(attr(xc, "filename"))
        if (grepl("uv", bf)) {
          x <- rwind::rotax(x, alfa = as.numeric(as.character(df$valfa[1])))
          if (grepl("pres", bf)) {
            x <- x[df$ilon[1], df$ilat[1], df$ilevel[1],,df$icomponent[1]]
          } else if (grepl("surf", bf)) {
            x <- x[df$ilon[1], df$ilat[1],,df$icomponent[1]]
          }
        } else if (grepl("hgt", bf)) {
          x <- x[df$ilon[1], df$ilat[1], df$ilevel[1],]
        } else if (grepl("mslp", bf)) {
          x <- x[df$ilon[1], df$ilat[1],] / 100
        } else if (grepl("omega", bf)) {
          x <- x[df$ilon[1], df$ilat[1], df$ilevel[1],]
        }
        cname <- if  (grepl("log", bf)) paste0("log_", varname) else varname
        if (grepl("uv", bf)) {
          if (grepl("surf", bf)) {
            cname <- paste0(cname, "_surf")
          }
        }
        if (is.null(ret)) {
          ret = matrix(x, ncol = 1)
          colnames(ret)[1] <- cname
          rownames(ret) <- names(x)
        } else {
          ret <- cbind(ret, x)
          colnames(ret)[ncol(ret)] <-  cname
        }

      }, error = function(e) {
        print(e)
      }, finally = {
        next
      })
    }
  }
  return(ret)
}

save_data_by_cor <- function(fun = max, file = "out.csv") {
  data <- get_data_by_cor(fun)
  save_csv(data, file)
}

cor2 <- function(x, pm, alfa = seq(0, 360, 20), par = T) {
  calc <- function(uv) {
    r <- if (is.null(dim(uv))) matrix(uv, ncol = 1) else uv
    rot <- function(a) rotax(uv, alfa = a)
    cor2 <- function(x1) apply(pm, 2, cor, y = x1, use = "pairwise.complete.obs")
    if (ncol(r) == 2) {
      r <- vapply(alfa, rot, outer(1:dim(r)[1], 1:dim(r)[2]))
      ret <- apply(r, 2:3, cor2)
      dimnames(ret)[[3]] <- alfa
    } else {
      ret <- apply(r, 2, cor2)
    }
    return(ret)
  }

  if (is.null(dim(pm))) pm <- matrix(pm, ncol = 1)
  dims <- (1:length(dim(x)))
  comp_loc <- which(dim(x) == 2) # find uv dim
  if (length(comp_loc) > 0) {
    dims <- dims[-comp_loc]
    dim_comp <- dim(x)[comp_loc]
    dim_comp_names <- dimnames(x)[[comp_loc]]
    dim_alfa <- length(alfa)
  } else {
    alfa = 0
    dim_comp <- 1
    dim_comp_names <- "X"
    dim_alfa <- 1
  }
  time_loc <- which(dim(x) == nrow(pm)) # find time dim
  if (length(time_loc) == 0)
    stop("x does not have a dim suitable with pm")
  dims <- dims[-time_loc]
  if (require(rwind)) {
    if (require(parallel) && par) {
      cl <- parallel::makeCluster(getOption("cl.cores", detectCores() / 2))
      parallel::clusterExport(cl, c("pm", "alfa"), envir = environment())
      parallel::clusterEvalQ(cl, library(rwind))
      r <- parallel::parApply(cl, x, dims, calc)
      parallel::stopCluster(cl)
    } else {
      r <- apply(x, dims, calc)
    }
  } else {
    stop("Install rwind package")
  }
  r <- array(r, dim = c(ncol(pm), dim_comp,
                        dim_alfa, dim(x)[dims]))
  names(dim(r)) <- c("pm",
                     "component",
                     "alfa",
                     names(dim(x)[dims]))
  dimnames(r) <- append(list(pm = colnames(pm),
                             component = dim_comp_names,
                             alfa = alfa), dimnames(x)[dims])
  dim_order <- 1:length(dim(r))
  first <- match(dim(x)[dims], dim(r))
  remain <- dim_order[-first]
  r <- aperm(r, c(first, remain))

  # copy attributes
  attrs <- names(attributes(x))
  attrs <- attrs[!(attrs %in% c("dim", "dimnames"))]
  for (a in attrs) attr(r, a) <- attr(x, a)

  return(drop(r))
}

