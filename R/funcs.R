#' Common functions
#' 
#' push an object to the list
#' 
#' @param l, that are going to push to the list
#' @param ..., list object
#' @return combined list
#' 
#' @export
#' 
#' @examples
#'    mylist <- list()
#'    newlist <- push ( 1, mylist )

push <- function(l, ...) c(l, list(...))
