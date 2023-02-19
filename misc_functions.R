# Fuction equivalent to dplyr's spread function, taking a long-form table and
# converting to a wide-form, with rows and columns as specified by rows and
# cols. If values is NULL, resulting table is populated by zeros and ones.
# If values is a numeric vector, the table is populated by the values aggregated
# by the supplied function (default sum)

ct <- function (rows, 
                cols, 
                values = NULL, 
                FUN = sum, 
                convertNAToZero = TRUE,...) 
{
  if(!is.vector(rows)) rows <- as.vector(rows)
  if(!is.vector(cols)) cols <- as.vector(cols)
  if(is.null(values)) values <- rep(1,length(rows))
  results <- tapply(values, list(rows, cols), FUN, ...)
  if(convertNAToZero)
    results[is.na(results)] <- 0
  results
}

unscale <- function(y_scaled_transformed, x, scaled = TRUE, log = TRUE, log_add = 0){  
  #y_scaled_transformed is modelled output of variable x calculated on scaled 
  # (if scale = TRUE), and logged (if log = TRUE) data, x
  # x is the original vector of raw modelled data before transformation (if log = TRUE) and scaling
  # log_add is the amount added before logging to avoid log(0)
  if(log){
    x_transformed <- log(x + log_add)
  }else{
    x_transformed <- x
  }
  if(scaled){
    scale_x <- scale(x_transformed) 
    y_transformed <- y_scaled_transformed * attr(scale_x, 'scaled:scale') + 
      attr(scale_x, 'scaled:center')
  }else{
    y_transformed <- y_scaled_transformed
  }
  if(log){
    y <- exp(y_transformed) - log_add
  }else{
    y <- y_transformed
  }
  list(y_transformed = y_transformed, 
       y = y)
}
