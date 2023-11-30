#' Akaike Information Criterion (AIC) of a \code{tvcure.object}.
#'
#' @description
#' Akaike Information Criterion (AIC) for the fitted tvcure model in a \code{tvcure.object}.
#'
#' @usage \method{AIC}{tvcure}(x, ..., k=2)
#'
#' @param x an object of class \code{\link{tvcure.object}}.
#' @param k the penalty per parameter to be used (Default: k=2 for the classical AIC).
#' @param ... optionally more fitted objects.
#'
#' @details Provides the AIC of the fitted tvcure model in a given \code{\link{tvcure.object}}.
#'
#' @return The AIC of the tvcure model in \code{x}.
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#' @references Lambert, P. and Kreyenfeld, M. (2023). Exogenous time-varying covariates in double additive cure survival model
#' with application to fertility.
#' \emph{Journal of the Royal Statistical Society, Series A}, in press.
#'
#' @examples
#' require(tvcure)
#' ## data(tvcure_Data)
#' ## fit = tvcure(...)
#' ## AIC(fit)
#'
#' @seealso \code{\link{tvcure}}, \code{\link{tvcure.object}}, \code{\link{BIC.tvcure}}
#'
#' @export
#'
AIC.tvcure <- function(x, ..., k=2){
    obj = x
    lls = function(obj) return(ans = c(dev=obj$fit$dev, edf=obj$fit$ED.tot, nobs=obj$fit$nobs))
    if (!missing(...)) {
        vals = sapply(list(obj,...), lls)
        val <- data.frame(edf = round(vals[2L, ],2), AIC = vals[1L, ] + k * vals[2L, ])
        nos <- na.omit(vals[3L, ])
        if (length(nos) && any(nos != nos[1L])) warning("models are not all fitted to the same number of observations")
        Call <- match.call()
        Call$k <- NULL
        row.names(val) <- as.character(Call[-1L])
        val
    } else {
        vals = unname(lls(obj))
        vals[1L] + k * vals[2L]
    }
}
