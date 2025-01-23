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

RegFxn <- function(x, y, q){
  if(frequency(x) == 12){
    midas_r(y ~ mls(y, 1, 1) + mls(x, 12 + 0:11, 12, nealmon),
            start = list(
              x = nealmon_start_params(y, x, 11, 12, q)
            )
    )
  } else {
    midas_r(y ~ mls(y, 1, 1) + mls(x, 4 + 0:3, 4, nealmon),
            start = list(
              x = nealmon_start_params(y, x, 3, 4, q)
            )
    )
  }
}

rmseWeights <- function(Accuracy){
  w <- matrix(data = NA, ncol = 1)
  for(i in 1:nrow(Accuracy)){
    w[i] <- ( 1 / ( Accuracy[i, 2] ^ 2 ) ) / sum( 1 / ( Accuracy[, 2] ^ 2 ) ) 
  }
  return(w)
}

Forecast <- function(x, y, q, h){
  y_train <- window(y, end = end(time(y))[1] - h)
  y_test <- y
  
  models <- lapply(x, RegFxn, y = y_train, q = q)
  
  fc <- mapply(
    function(model, x) #NewData
      forecast(
        model,
        se = TRUE,
        level = c(68, 95),
        method = 'static',
        #npaths = 3000,
        insample = get_estimation_sample(model),
        newdata = list(
          x = rep(NA, (frequency(x) * h)) #rnorm(frequency(NewData) * h, mean(NewData), sd(NewData))
        )
      ),
    models,
    x,
    SIMPLIFY = FALSE
  )
  
  Accuracy <- data.frame(
    variable = names(x),
    'RMSE' = sapply(fc, function(z) accuracy(z, y_test)[2, 'RMSE']),
    row.names = 1:length(x)
  )
  
  W <- rmseWeights(Accuracy)
  
  Forecast_STG <- data.frame(
    'mean' = sapply(fc, function(a) cbind(summary(a)$mean)),
    'Lower95' = sapply(fc, function(a) cbind(summary(a)$lower))[2, ],
    'Lower68' = sapply(fc, function(a) cbind(summary(a)$lower))[1, ],
    'Upper68' = sapply(fc, function(a) cbind(summary(a)$upper))[1, ],
    'Upper95' = sapply(fc, function(a) cbind(summary(a)$upper))[2, ]
  )
  
  weighted_fc <- t(sapply(Forecast_STG, function(b) b %*% W)) %>% ts(., end = end(time(y))[1])
  
  out <- structure(
    list(
      mean = weighted_fc[, 'mean'],
      lower = cbind(
        lower68 = weighted_fc[, 'Lower68'],
        lower95 = weighted_fc[, 'Lower95']
      ),
      upper = cbind(
        upper68 = weighted_fc[, 'Upper68'],
        upper95 = weighted_fc[, 'Upper95']
      ),
      x = y_train,
      level = c(68, 95),
      method = 'midas_r'
    ),
    class = 'forecast'
  )
  return(out)
}