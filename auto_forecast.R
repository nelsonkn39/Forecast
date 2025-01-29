auto_forecast <- function(x, y, q, h, yk, xk, na = c(TRUE, FALSE), error = c("ME", "RMSE", "MAE", "MPE", "MAPE"), se = c(TRUE, FALSE), level = c(95, 80), method = c('static', 'dynamic'), npaths, newdata = NULL){
  
  y_freq <- frequency(y)
  if (y_freq == 1) {
    y_train <- window(y, end = end(y)[1] - h)
  } else if (h %% y_freq == 0) {
    y_train <- window(y, end = c(end(y)[1] - (y_freq / h), end(y)[2]))
  } else if (h > y_freq) {
    y_train <- window(y,  end = c((end(y)[1] - (h %/% y_freq)), (end(y)[2] - (h %% y_freq))))
  } else if (h < y_freq) {
    if ((h %% y_freq) >= end(y)[2]) {
      y_train <- window(y, end = c(end(y)[1] - 1, end(y)[2] + 1))
    } else {
      y_train <- window(y, end = c(end(y)[1], end(y)[2] - h))
    }
  }
  y_test <- y
  
  models <- lapply(x, reg_fn, y = y_train, q = q, yk = yk, xk = xk, na = na, h = h)
  
  if(na == FALSE){
    fc <- mapply(
      function(model, NewData)
        forecast(
          model,
          se = se,
          level = level,
          method = method,
          npaths = npaths,
          insample = get_estimation_sample(model),
          newdata = list(
            x = rnorm(frequency(NewData) / (frequency(NewData) / frequency(y)) * h, mean(NewData), sd(NewData))
          )
        ),
      models,
      x,
      SIMPLIFY = FALSE
    )
  } else {
    fc <- mapply(
      function(model, NewData)
        forecast(
          model,
          se = se,
          level = level,
          method = method,
          npaths = npaths,
          insample = get_estimation_sample(model),
          newdata = list(
            x = rep(NA, (frequency(NewData) / frequency(y)) * h)
          )
        ),
      models,
      x,
      SIMPLIFY = FALSE
    )
  }
  
  Accuracy <- data.frame(
    variable = names(x),
    error = sapply(fc, function(z) accuracy(z, y)[2, error]),
    row.names = 1:length(x)
  )
  
  W <- forecast_weights(Accuracy)
  
  forecast_df <- data.frame(
    sapply(fc, function(a) cbind(summary(a)$mean)) %*% W,
    sapply(fc, function(a) cbind(summary(a)$lower)[, 1]) %*% W,
    sapply(fc, function(a) cbind(summary(a)$lower)[, 2]) %*% W,
    sapply(fc, function(a) cbind(summary(a)$upper)[, 2]) %*% W,
    sapply(fc, function(a) cbind(summary(a)$upper)[, 1]) %*% W
  )
  colnames(forecast_df) <- c('Point Forecast', paste(c(rep('Lower', 2), rep('Upper', 2)), rep(level, 2)))
  
  if (y_freq == 1) {
    forecast_df <- ts(forecast_df, end = end(y))
  } else {
    forecast_df <-  ts(forecast_df, frequency = frequency(y), end = c(end(y)[1], end(y)[2]))
  }
  
  out <- structure(
    list(
      models = models,
      forecasts = fc,
      mean = forecast_df[, 1],
      lower = cbind(
        lower_a = forecast_df[, 2],
        lower_b = forecast_df[, 3]
      ),
      upper = cbind(
        upper_a = forecast_df[, 4],
        upper_b = forecast_df[, 5]
      ),
      x = y_train,
      level = level,
      method = 'midas_r'
    ),
    class = 'forecast'
  )
  return(out)
}