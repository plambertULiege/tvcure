#' Bayesian Information Criterion (BIC) of a \code{tvcure.object}.
#'
#' @description
#' Bayesian Information Criterion (BIC) for the fitted tvcure model in a \code{tvcure.object}.
#'
#' @usage \method{BIC}{tvcure}(x, ...)
#'
#' @param x An object of class \code{\link{tvcure.object}}.
#' @param ... Optionally more fitted objects.
#'
#' @details Bayesian (Schwarz) information criterion in a tvcure object, with a penalty calculated using the total effective degrees of freedom and the total number of observed events, -2log(L) + log(d)*ED.tot, smaller values being preferred during model selection.
#'
#' @return The BIC of the fitted tvcure model in \code{x}.
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#' @references Lambert, P. and Kreyenfeld, M. (2024). Exogenous time-varying covariates in double additive cure survival model
#' with application to fertility. \emph{Journal of the Royal Statistical Society, Series A}, in press.
#'
#' @examples
#' require(tvcure)
#' ## data(tvcure_Data)
#' ## fit = tvcure(...)
#' ## BIC(fit)
#'
#' @seealso \code{\link{tvcure}}, \code{\link{tvcure.object}}, \code{\link{AIC.tvcure}}
#'
#' @export
#' 
BIC.tvcure <- function(x, ...){
    obj = x
    lls = function(obj) return(ans = c(dev=obj$fit$dev, edf=obj$fit$ED.tot, d=obj$fit$d))
    ## lls = function(obj) return(ans = c(dev=obj$fit$dev, edf=obj$fit$ED.tot, nobs=obj$fit$nobs))
    if (!missing(...)) {
        vals = sapply(list(obj,...), lls)
        val <- data.frame(edf = round(vals[2L, ],2), BIC = vals[1L, ] + log(vals[3L, ]) * vals[2L, ])
        nos <- na.omit(vals[3L, ])
        if (length(nos) && any(nos != nos[1L])) warning("models are not all fitted to the same number of observations")
        Call <- match.call()
        row.names(val) <- as.character(Call[-1L])
        val
    } else {
        vals = unname(lls(obj))
        vals[1L] + log(vals[3L]) * vals[2L]
    }
}
