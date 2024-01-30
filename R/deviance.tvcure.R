#' Deviance of a \code{tvcure.object}.
#'
#' @description
#' Deviance for the fitted tvcure model in a \code{tvcure.object}.
#'
#' @usage \method{deviance}{tvcure}(object, ...)
#'
#' @param object An object of class \code{\link{tvcure.object}}.
#' @param ... Optionally more fitted objects.
#'
#' @details Deviance -2log(L/L.hat)
#'
#' @return The deviance of the fitted tvcure model in \code{x}.
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#' @references Lambert, P. and Kreyenfeld, M. (2024). Exogenous time-varying covariates in double additive cure survival model
#' with application to fertility. \emph{Journal of the Royal Statistical Society, Series A}, under review.
#'
#' @examples
#' require(tvcure)
#' ## Simulated data generation
#' beta = c(beta0=.4, beta1=-.2, beta2=.15) ; gam = c(gam1=.2, gam2=.2)
#' data = simulateTVcureData(n=500, seed=123, beta=beta, gam=gam,
#'                           RC.dist="exponential",mu.cens=550)$rawdata
#' ## TVcure model fitting
#' tau.0 = 2.7 ; lambda1.0 = c(40,15) ; lambda2.0 = c(25,70) ## Optional
#' model = tvcure(~z1+z2+s(x1)+s(x2), ~z3+z4+s(x3)+s(x4), data=data,
#'                tau.0=tau.0, lambda1.0=lambda1.0, lambda2.0=lambda2.0)
#' deviance(model)
#'
#' @seealso \code{\link{tvcure}}, \code{\link{tvcure.object}}, \code{\link{AIC.tvcure}}, \code{\link{BIC.tvcure}}, \code{\link{levidence}}
#'
#' @export
#'
deviance.tvcure <- function(object, ...){
    obj = object
    lls = function(obj) return(ans = c(dev=obj$fit$dev, edf=obj$fit$ED.tot, d=obj$fit$d))
    if (!missing(...)) {
        vals = sapply(list(obj,...), lls)
        val <- data.frame(edf = round(vals[2L, ],2), dev = vals[1L, ])
        nos <- na.omit(vals[3L, ])
        if (length(nos) && any(nos != nos[1L])) warning("models are not all fitted to the same number of observations")
        Call <- match.call()
        row.names(val) <- as.character(Call[-1L])
        val
    } else {
        vals = unname(lls(obj))
        vals[1L]
    }
}
