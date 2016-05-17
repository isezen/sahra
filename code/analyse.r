# Saharan Dust Transport Research
# 2016-05-16 Ismail SEZEN
# sezenismail@gmail.com

source("code/correlation.r")
source("code/filehelper.r")

reg_tree <- function(x, ...) {
  d <- x
  if (!is.formula(d)) {
    stopifnot(!is.null(dim(d)) & length(dim(d)) == 2)
    d <- as.data.frame(d)
  }
  fit <- rpart(d, na.action = na.rpart, control = rpart.control(minsplit = 20))
  par(xpd = TRUE)
  plot(fit, uniform = TRUE, col = "gray", ...)
  text(fit, use.n = T, all = T, pos = 1,
       fancy = T, fheight = 4, fwidth = 7, bg = "antiquewhite")
  return(fit)
}

regtree_analyse <- function() {
  pm <- as.matrix(read_pm10()[, 6]) # TotalAveragePM_24h
  hgt <- as.matrix(read_hgt()[, 1]) # HGHT1
  data <- get_data_by_cor()[,6:10] # do not get log results
  x <- cbind(pm, data[,c(1,2,3,4,5), drop = F])
  pdf("reg_tree.pdf", width = 16, height = 16)
  fit <- reg_tree(x)
  dev.off()
  leafs <- xts(fit$where, order.by = as.POSIXct(names(fit$where), tz = "GMT"))
  leafs <- merge(leafs, zoo(, as.POSIXct(rownames(pm), tz = "GMT")), all = T)
  x2 <- cbind(x, as.matrix(leafs))
  # names(fit$where[fit$where == split_id[1]])
  cols <- rep(palette(), 10)
  print(cols[unique(leafs)])

  open3d()
  cols <- c("black", "green3", "blue", "cyan", "magenta", "red")
  plot3d(x[,c(5,6,1)], col = cols[as.numeric(as.factor(leafs))],
         type = "s", size = 0.5)

  open3d()
  cols <- c(rep("gray90", 4), "red", "red")
  plot3d(x[,c(5,6,1)], col = cols[as.numeric(as.factor(leafs))],
         type = "s", size = 0.5)
  rgl.bbox(xlen = 0, ylen = 0, zlen = 0, color = c('gray50'))
  varname <- "HGHT1"

  plot(y = x2[,1], x = x2[,varname], pch = 21, las = 1,
       ylab = "log_pm10", xlab = varname)
  for (p in unique(leafs)) {
    if (!is.na(p)) {
      df <- x2[x2[,"leafs"] == p,]
      points(y = df[,1], x = df[,varname], pch = 21, bg = p)
    }
  }
}
