#' Test package availability
#'
#' @description Question the user on core International Geo Sample Number, sampling number, etc.
#'
#' @export
#' @examples pkgTest('vegan')
#'

pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}
