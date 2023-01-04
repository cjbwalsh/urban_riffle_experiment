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