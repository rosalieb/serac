#' Round values
#'
#' @description Round value at the upper or lower value, defined by 'base'. Used for plotting purpose.
#'
#' @export
#' @param x number to round up
#' @param base round up exponant
#' @keywords visualisation
#' @examples mround(142,10)
#'

mround <- function(x, base)
{
  if (x>0) base*ceiling(x/base) else base*floor(x/base)
}
