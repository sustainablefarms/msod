#' @title Converting BUGS Variable Names
#' @description For converting values for array-valued parameters from the bugs variable format to an array
#' @param values is a list of values named according to the bugs variables names
#' @param varname is the desired variable name (e.g. 'u.b')
#' @param rowidx is a list of rows to extract, by number
#' @param colidx is a list of columns to extract, by number
#' @export
bugsvar2array <- function(values, varname, rowidx, colidx){
  if (is.vector(values)) {
    values <- matrix(values, nrow = 1, dimnames = list(row = NULL, col = names(values)))
  }
  idx <- expand.grid(row = rowidx, col = colidx)
  bugsnames <- paste0(varname, "[",idx$row, ",", idx$col, "]") #order matters, expand.grid must go through rows and then columns
  value <- array(t(values[, bugsnames]), 
                 dim = c(length(rowidx), length(colidx), nrow(values)), 
                 dimnames = list(row = rowidx, col = colidx, draw = 1:nrow(values)))
  return(value)
}

#' @describeIn bugsvar2array For converting values for array-valued parameters from the bugs variable format to a matrix
#' @param values is a list of values named according to the bugs variables names
#' @param varname is the desired variable name (e.g. 'u.b')
#' @param rowidx is a list of rows to extract, by number
#' @param colidx is a list of columns to extract, by number
#' @export
bugsvar2matrix <- function(values, varname, rowidx, colidx){
  arr <- bugsvar2array(values, varname, rowidx, colidx)
  mat <- drop_to_matrix(arr, 3)
  return(return(mat))
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

#' @title A helper function to get a vector of parameters from a fitted object
#' @param fit fitted runjags object with summary included
#' @param type An integer will select the draw from the posterior
#' (e.g. type = 5 will return the vector of parameters in the 5th sample from the posterior, mcmc chains are concatenated)
#' A charactor vector can be used to the parameters based on summarary statistics like quantiles and moments.
#' Supported values so far "median" and "mean".
#' If [type] is "marginal" or "all" then all draws are returned.
#' @export
get_theta <- function(fit, type){
  if (is.numeric(type) && length(type) > 1 && !is.null(names(type))){ #assumed passed 'type' is actually the desired theta
    return(type)
  }
  if (is.numeric(type) && length(type) == 1){
    chainidx <- floor(type / (fit$sample + 1 )) + 1
    sampleinchain <- type - fit$sample * (chainidx - 1)
    theta <- fit$mcmc[[chainidx]][sampleinchain, ]
    return(theta)
  }
  if (type == "median"){theta <- apply(do.call(rbind, fit$mcmc), MARGIN = 2, median)}
  if (type == "mean"){theta <- apply(do.call(rbind, fit$mcmc), MARGIN = 2, mean)}
  if (type == "marginal" | type == "all"){theta <- do.call(rbind, fit$mcmc)}
  return(theta)
}

#' @title A quick replacement to [runjags::list.format()] that does nothing if the argument is already a list.
#' @param data Same as [runjags::list.format()]. Typically found in the `data` slot of runjags object.
#' @param checkvalid See [runjags::list.format()].
#' @export as_list_format
as_list_format <- function(data, checkvalid = TRUE){
  if ("list" %in% class(data)){return(data)}
  out <- runjags::list.format(data, checkvalid = checkvalid)
  return(out)
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