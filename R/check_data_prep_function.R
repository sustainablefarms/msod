#' @examples
#' myfun <- function(B, sd){
#'   out <- (B - cmns) * sd
#'   return(out)
#' }
#' runsusingpackagesonly(myfun, args = list(sd = 3)) #returns FALSE
#' B <- matrix(rnorm(4 * 3), nrow = 3)
#' runsusingpackagesonly(myfun, args = list(sd = 3, B = B)) #returns FALSE
#' cmns <- colMeans(B) 
#' runsusingpackagesonly(myfun, args = list(sd = 3, B = B)) #returns FALSE too - myfun run such that it can't see globalenv()
#'     
#' myfun2 <- function(B, sd, cmns){
#'   out <- (B - cmns) * sd
#'   return(out)
#' }
#' runsusingpackagesonly(myfun2, args = list(sd = 3, B = B, cmns = cmns)) # TRUE, which is correct! :)

#' @title Check that function runs using objects only in its argument list or imported from currently loaded packages
#' @param fun The function to check
#' @param args An argument list to pass to `fun` like in [`do.call()`].
#' @return TRUE if the functions runs in the parent environment of the global environment (an environment with only packages) with arguments attached, 
#' which means the function runs when it can't access the interactively defined objects.
#' Or FALSE if arguments are missing, the function tries to access an object in the global environment (and not in the parent environment), 
#' or if there is some other problem with the function executing. In this case a warning is given. This warning contains the error condition that created the failure.
#' @export
runsusingpackagesonly <- function(fun, args, erroronfalse = TRUE){
  stop("Function obsolete")
  environment(fun) <- parent.env(globalenv())  #an environment with all packages loaded, but none of the interactively defined objects
  out <- try(do.call(fun, args), silent = TRUE) #run the function
  if (!("try-error" %in% class(out))){return(TRUE)}
  if (erroronfalse){
    stop(attr(out, "condition"))
  } else {
    warning(attr(out, "condition")) 
    return(FALSE)
  }
}

