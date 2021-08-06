#' @title Converting Parameter Values in BUGS Name and Format to Arrays and Matrices
#' @description For converting values for array-valued parameters from the bugs variable format to an array
#' @param values is a list of values named according to the bugs variables names
#' @param varname is the desired variable name (e.g. 'occ.b')
#' @param rowidx is a list of rows to extract, by number
#' @param colidx is a list of columns to extract, by number
#' @export
bugsvar2array <- function(values, varname, rowidx, colidx){
  if (is.vector(values)) {
    values <- matrix(values, nrow = 1, dimnames = list(row = NULL, col = names(values)))
  }
  #checks
  nvalues <- sum(grepl(paste0("^", varname, "\\[.*"), colnames(values)))
  stopifnot(nvalues > 0)
  stopifnot(length(rowidx) * length(colidx) <= nvalues)
  #actual extraction
  idx <- expand.grid(row = rowidx, col = colidx)
  bugsnames <- paste0(varname, "[",idx$row, ",", idx$col, "]") #order matters, expand.grid must go through rows and then columns
  value <- array(t(values[, bugsnames]), 
                 dim = c(length(rowidx), length(colidx), nrow(values)), 
                 dimnames = list(row = rowidx, col = colidx, draw = 1:nrow(values)))
  return(value)
}

#' @describeIn bugsvar2array For converting values for array-valued parameters from the bugs variable format to a matrix
#' @param values is a list of values named according to the bugs variables names
#' @param varname is the desired variable name (e.g. 'occ.b')
#' @param rowidx is a list of rows to extract, by number
#' @param colidx is a list of columns to extract, by number
#' @export
bugsvar2matrix <- function(values, varname, rowidx, colidx){
  arr <- bugsvar2array(values, varname, rowidx, colidx)
  mat <- drop_to_matrix(arr, 3)
  return(return(mat))
}

#' @describeIn bugsvar2array Converting vector-valued parameters from the bugs variable format to an R array
#' @param idx Index for the vector
#' @export
bugsvar2array_vector <- function(values, varname, idx){
  if (is.vector(values)) {
    values <- matrix(values, nrow = 1, dimnames = list(row = NULL, col = names(values)))
  }
  #checks
  nvalues <- sum(grepl(paste0("^", varname, "\\[.*"), colnames(values)))
  stopifnot(nvalues > 0)
  stopifnot(length(idx) <= nvalues)
  #actual extraction
  bugsnames <- paste0(varname, "[",idx, "]")
  value <- array(t(values[, bugsnames]),
        dim = c(length(idx), nrow(values)),
        dimnames = list(idx = idx, draw = 1:nrow(values)))
  return(value)  
}

#' @describeIn bugsvar2array Converts a matrix of parameter values to bugs variable format
#' @param theta is a parameter arrays
#' @param name parameter
#' @return a named vector. Names are given by name and the index in the array
#' @export
matrix2bugsvar <- function(theta, name){
  values <- as.vector(theta) #runs through first dimension, then second dimension, then third dimension...
  idx <- expand.grid(row = 1:nrow(theta), col = 1:ncol(theta))
  bugsnames <- paste0(name, "[",idx$row, ",", idx$col, "]") #order matters, expand.grid must go through rows and then columns
  names(values) <- bugsnames
  return(values)
}

drop_to_matrix <- function(array, dimdrop = 3){
  stopifnot(dim(array)[dimdrop] == 1) #must use subsetting [,,1, drop = FALSE] first
  if ((dim(array))[[1]] == 1){
    return(matrix(array, nrow = 1))
  } else if ((dim(array))[[2]] == 1){
    return(matrix(array, ncol = 1))
  } else {
    return(drop(array))
  }
}