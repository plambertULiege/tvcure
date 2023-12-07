#' Log-evidence of a tvcure object.
#'
#' @description
#' The log-evidence of the fitted tvcure model in a \code{tvcure.object}.
#'
#' @usage levidence(x, ...)
#'
#' @param x A \code{\link{tvcure.object}}.
#' @param ... Optionally more tvcure objects.
#'
#' @details Provides the log-evidence (or log-marginal likelihood) of the fitted tvcure model in a given \code{\link{tvcure.object}}, where the evidence is the marginal posterior of the penalty parameters at their selected values.
#' @return The log-evidence of the tvcure model in \code{x}.
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#' @references Lambert, P. and Kreyenfeld, M. (2024). Exogenous time-varying covariates in double additive cure survival model
#' with application to fertility. \emph{Journal of the Royal Statistical Society, Series A}, under review.
#'
#' @examples
#' require(tvcure)
#' ## data(tvcure_Data)
#' ## fit = tvcure(...)
#' ## levidence(fit)
#'
#' @seealso \code{\link{tvcure}}, \code{\link{tvcure.object}}, \code{\link{AIC.tvcure}}, \code{\link{BIC.tvcure}}
#'
#' @export
#'
#' @examples
#' require(tvcure)
#' ## Simulated data generation
#' beta = c(beta0=.4, beta1=-.2, beta2=.15) ; gam = c(gam1=.2, gam2=.2) 
#' df.raw = simulateTVcureData(n=500, seed=123, beta=beta, gam=gam,
#'                           RC.dist="exponential",mu.cens=550)$df.raw
#' ## TVcure model fitting
#' tau.0 = 2.5 ; lambda1.0 = c(285,15) ; lambda2.0 = c(25,1325) ## Optional
#' model = tvcure(~z1+z2+s(x1)+s(x2), ~z3+z4+s(x3)+s(x4), df=df.raw,
#'                tau.0=tau.0, lambda1.0=lambda1.0, lambda2.0=lambda2.0)
#' levidence(model)

##levidence <- function(x, ...) UseMethod("levidence")
## levidence.tvcure <- function(x, ...){
levidence <- function(x, ...){
    obj = x
  lls = function(obj) return(ans = c(levidence=obj$fit$levidence, edf=obj$fit$ED.tot, nobs=obj$fit$nobs))
  if (!missing(...)) {
    vals = sapply(list(obj,...), lls)
    val <- data.frame(edf = round(vals[2L, ],2), levidence = vals[1L, ])
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