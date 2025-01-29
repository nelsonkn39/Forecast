reg_fn <- function(x, y, q, h, yk, xk, na = c(TRUE, FALSE)){
  y_freq <- frequency(y)
  x_freq <- frequency(x)
  if (na == FALSE) {
    model <- midas_r(y ~ mls(y, 1:yk, 1) + fmls(x, xk, (x_freq / y_freq), nealmon),
                     start = list(
                       x = nealmon_start_params(y, x, xk, (x_freq / y_freq), q)
                     )
    )
  } else {
    model <- midas_r(y ~ mls(y, 1:yk, 1) + mls(x, (((x_freq / y_freq) * h) + 0:xk), (x_freq / y_freq), nealmon),
                     start = list(
                       x = nealmon_start_params(y, x, xk, (x_freq / y_freq), q)
                     )
    )
  }
  return(model)
}