recode <- function(X, items = NULL, values = defaultValues){
  defaultValues = min(X, na.rm = TRUE) : max(X, na.rm = TRUE)
  X[, items][!is.na(X[, items])] <- max(values) - X[, items][!is.na(X[, items])] + min(values)
  return(X)
}
