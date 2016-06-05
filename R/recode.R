recode <- function(X, items = NULL, values = defaultValues){
  defaultValues = min(X) : max(X)
  X[, items][!is.na(X[, items])] <- max(values) - X[, items][!is.na(X[, items])] + min(values)
  return(X)
}
