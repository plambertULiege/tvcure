#' Specification of smooth terms in formulas in the tvcure function.
#' @keywords internal
#' @export
#' @param x Name of the variable for which an additive term is requested.
#' @param (Optional) reference value for \code{x} where the additive term is zero.
#' @return The submitted variable for which an additive term is required.
s <- function(x,ref=NULL) x
