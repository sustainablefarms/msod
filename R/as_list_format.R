#' @title A quick replacement to [runjags::list.format()] that does nothing if the argument is already a list.
#' @param data Same as [runjags::list.format()]. Typically found in the `data` slot of runjags object.
#' @param checkvalid See [runjags::list.format()].
#' @export as_list_format
as_list_format <- function(data, checkvalid = TRUE){
  if ("list" %in% class(data)){return(data)}
  out <- runjags::list.format(data, checkvalid = checkvalid)
  return(out)
}
