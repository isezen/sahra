# Saharan Dust Transport Research
# 2016-05-04 Ismail SEZEN
# sezenismail@gmail.com

# http://stackoverflow.com/questions/14500707/select-along-one-of-n-dimensions-in-array
index_array <- function(x, dim, value, drop = T) {
  # Create list representing arguments supplied to [
  # bquote() creates an object corresponding to a missing argument
  indices <- rep(list(bquote()), length(dim(x)))
  indices[[dim]] <- value

  # Generate the call to [
  call <- as.call(c(
    list(as.name("["), quote(x)),
    indices,
    list(drop = drop)))
  # Print it, just to make it easier to see what's going on
  # print(call)

  # Finally, evaluate it
  eval(call)
}
