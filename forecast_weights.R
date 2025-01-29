forecast_weights <- function(Accuracy){
  w <- matrix(data = NA, ncol = 1)
  for(i in 1:nrow(Accuracy)){
    w[i] <- ( 1 / ( Accuracy[i, 2] ^ 2 ) ) / sum( 1 / ( Accuracy[, 2] ^ 2 ) ) 
  }
  return(w)
}