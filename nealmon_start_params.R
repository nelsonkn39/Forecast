nealmon_start_params <- function(y, x, k, m, q) {
  d <- k + 1
  i <- 1:d
  p <- poly(i, degree = q - 1, raw = TRUE)
  
  mt <- matrix(
    fmls(x, k, m) %*% p[, 1:ncol(p)],
    ncol = q - 1,
    dimnames = list(1:length(y), paste0('z', 1:(q - 1)))
  )
  
  z_s <- ts(
    mt,
    start = start(y)
  )
  
  zmod <- lm(y ~ z_s)
  
  return(as.vector(coef(zmod)))
}