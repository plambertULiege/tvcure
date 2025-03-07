#' Generic function for computing log-evidence
#'
#' This is a generic function returning the log-evidence for the class of the input object.
#'
#' @param object An object for which log-evidence is to be computed.
#' @param ... Additional arguments (not used for now).
#' @export
logEvid <- function(object, ...) {
  UseMethod("logEvid")
}

#' Log-evidence of a tvcure object.
#'
#' @description
#' The log-evidence of the fitted tvcure model in a \code{tvcure.object}.
#'
#' @param object A \code{\link{tvcure.object}}.
#' @param ... Optionally more tvcure objects.
#'
#' @details Provides the log-evidence (or log-marginal likelihood) of the fitted tvcure model in a given \code{\link{tvcure.object}}, where the evidence is the marginal posterior of the penalty parameters at their selected values.
#' @return The log-evidence of the tvcure model in \code{object}.
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#' @references Lambert, P. and Kreyenfeld, M. (2025).
#' Time-varying exogenous covariates with frequently changing values in double additive cure survival model: an application to fertility.
#' \emph{Journal of the Royal Statistical Society, Series A}. <doi:10.1093/jrsssa/qnaf035>
#'
#' @examples
#' require(tvcure)
#' ## data(tvcure_Data)
#' ## fit = tvcure(...)
#' ## logEvid(fit)
#'
#' @seealso \code{\link{tvcure}}, \code{\link{tvcure.object}}, \code{\link{AIC.tvcure}}, \code{\link{BIC.tvcure}}
#'
#' @export
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
#' logEvid(model)

##logEvid <- function(x, ...) UseMethod("logEvid")
## logEvid.tvcure <- function(x, ...){
logEvid.tvcure <- function(object, ...){
    obj = object
  lls = function(obj) return(ans = c(logEvid=obj$fit$logEvid, edf=obj$fit$ED.tot, nobs=obj$fit$nobs))
  if (!missing(...)) {
    vals = sapply(list(obj,...), lls)
    val <- data.frame(edf = round(vals[2L, ],2), logEvid = vals[1L, ])
    nos <- na.omit(vals[3L, ])
    if (length(nos) && any(nos != nos[1L])) warning("models are not all fitted to the same number of observations")
    Call <- match.call()
    Call$k <- NULL
    row.names(val) <- as.character(Call[-1L])
    val
  } else {
    vals = unname(lls(obj))
    vals[1L]
  }
}
