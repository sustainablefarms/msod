#' @title Save and load a data preparation function
#' @description Save function without saving all of the global environment
#' @param fun A data preparation function with first argument `indf` and other fixed parameters.
#' @param checkwith is a an object used only for checking that the function does not depend on items in the global environment
#' @param params A named list of any parameters that specify the transformation that fun performs.
#' @export
save_process <- function(fun, checkwith, params = NULL){
  if (names(formals(fun))[[1]] != "indf"){stop("Name of first argument must be 'indf'")}
  if (!setequal(names(formals(fun))[-1], names(params))){stop("The function has different parameter names to the parameters supplied.")}
  
  environment(fun) <- parent.env(globalenv())  #an environment with all packages loaded, but none of the interactively defined objects
  out <- try(do.call(fun, args = c(list(indf = checkwith), params)), silent = TRUE) #run the function
  if ("try-error" %in% class(out)){
    stop(paste("Could not run in packages-only envirnoment: ", attr(out, "condition")))
  }
  
  # following checks that te result doesn't change each time the function is called (repeatable, a little bit)
  out2 <- try(do.call(fun, args = c(list(indf = checkwith), params)), silent = TRUE)
  stopifnot(identical(out, out2))
  
  return(list(
    fun = fun,
    params = params
  ))
}

#' @export
apply_saved_process <- function(saved, indf, envir = parent.env(globalenv())){
  fun <- saved$fun
  environment(fun) <- envir
  processed <- do.call(fun, c(list(indf = indf), saved$params), quote = TRUE)
  return(processed)
}

#' @export
apply_process <- function(fun, indf, params = NULL, envir = parent.env(globalenv())){
  environment(fun) <- envir
  processed <- do.call(fun, c(list(indf = indf), params), quote = TRUE)
  return(processed)
}
