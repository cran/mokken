"check.data" <-
function(X){
   if (data.class(X) != "matrix" && data.class(X) != "data.frame")
     stop("Data are not matrix or data.frame")
    matrix.X <- as.matrix(X)
    if (is.na(any(X))) stop("Missing values are not allowed")
    if (any(mode(matrix.X)!="numeric")) stop("Data must be numeric")
    if (any(matrix.X) < 0) stop("Data should be positive")
    matrix.X <- matrix.X - min(matrix.X)
    return(matrix.X)
}
