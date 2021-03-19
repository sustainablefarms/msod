colprods_arr <- function(arr){
  arr2 <- aperm(arr, c(1, 3, 2))
  dim(arr2) <- c(dim(arr)[[1]] * dim(arr)[[3]], dim(arr)[[2]])
  colsum2 <- Rfast::rowprods(arr2)
  dim(colsum2) <- dim(arr)[c(1, 3)]
  colsum <- colsum2
  dimnames(colsum) <- list(dimnames(arr)[[1]], dimnames(arr)[[3]])
}

rowprods_arr <- function(arr){
  arr2 <- arr #aperm(arr, c(1, 3, 2))
  dim(arr2) <- c(dim(arr)[[1]],  dim(arr)[[2]] * dim(arr)[[3]])
  colsum2 <- Rfast::colprods(arr2)
  dim(colsum2) <- dim(arr)[c(2, 3)]
  colsum <- colsum2
  dimnames(colsum) <- list(dimnames(arr)[[2]], dimnames(arr)[[3]])
  return(colsum)
}

