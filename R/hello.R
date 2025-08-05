#' Simple greeting function
#'
#' @param name The name to be greeted by. If not set, it greets the whole world.
#'
#' @returns Message
#' @export
#'
#' @examples hello("Nigel")
hello <- function(name = "world") {
  print( paste("Hello", name) )
}
