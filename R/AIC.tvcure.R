#' Akaike Information Criterion (AIC) of a tvcure object.
#'
#' @description
#' Akaike Information Criterion (AIC) for the fitted tvcure model in a \code{tvcure.object}.
#'
#' @usage \method{AIC}{tvcure}(object, ..., k=2)
#'
#' @param object A \code{\link{tvcure.object}}.
#' @param k The penalty per parameter to be used. (Default: k=2 for the classical AIC).
#' @param ... Other optional tvcure objects.
#'
#' @details Akaike information criterion for the fitted model in a tvcure object, with a penalty calculated using the total effective degrees of freedom, -2log(L) + 2*ED.tot, smaller values being preferred during model selection.
#'
#' @return The AIC as a numeric value, computed according to the model specified in the input \code{object}.
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#' @references Time-varying exogenous covariates with frequently changing values in double additive cure survival model: an application to fertility.
#' \emph{Journal of the Royal Statistical Society, Series A}. <doi:10.1093/jrsssa/qnaf035>
#'
#' @examples
#' \donttest{
#' require(tvcure)
#' ## Simulated data generation
#' beta = c(beta0=.4, beta1=-.2, beta2=.15) ; gam = c(gam1=.2, gam2=.2)
#' data = simulateTVcureData(n=500, seed=123, beta=beta, gam=gam,
#'                           RC.dist="exponential",mu.cens=550)$rawdata
#' ## TVcure model fitting
#' tau.0 = 2.7 ; lambda1.0 = c(40,15) ; lambda2.0 = c(25,70) ## Optional
#' model = tvcure(~z1+z2+s(x1)+s(x2), ~z3+z4+s(x3)+s(x4), data=data,
#'                tau.0=tau.0, lambda1.0=lambda1.0, lambda2.0=lambda2.0)
#' AIC(model)
#' }
#'
#' @seealso \code{\link{tvcure}}, \code{\link{tvcure.object}}, \code{\link{BIC.tvcure}}, \code{\link{logEvid}}
#'
#' @export
#'
AIC.tvcure <- function(object, ..., k=2){
    obj = object
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
