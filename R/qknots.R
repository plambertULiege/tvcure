#' Specification of the knots of a cubic B-spline basis for given data.
#'
#' @usage qknots(x, xmin=NULL, xmax=NULL,
#'        equid.knots = TRUE, pen.order=2, K=25)
#' @param x data that should be supported by the knots of the B-spline basis.
#' @param xmin (Optional) minimum value for the knots.
#' @param xmax (Optional) maximum value for the knots.
#' @param equid.knots Logical indicating if equidistant knots are desired (Default: TRUE).
#' @param pen.order penalty order (if equid.knots = TRUE) (Default: 2).
#' @param K number of B-splines in the basis (Default: 25).
#'
#' @return a list containing the following elements:
#' \itemize{
#' \item \code{xmin} : minimum value of the knots.
#' \item \code{xmax} : maximum value of the knots.
#' \item \code{knots} : vector containing the knots: equidistant if \code{equid.knots} is TRUE, based on quantiles of \code{x} otherwise.
#' \item \code{Pd} : penalty matrix for the B-spline coefficients.
#' \item \code{pen.order} : penalty order for the P-spline model.
#' }
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#' @references Lambert, P. and Kreyenfeld, M. (2024). Exogenous time-varying covariates in double additive cure survival model
#' with application to fertility. \emph{Journal of the Royal Statistical Society, Series A}, under review.
#'
#' @export
#'
#' @examples
#' x = rnorm(100)
#' qknots(x)
qknots = function(x, xmin=NULL,xmax=NULL, equid.knots = TRUE, pen.order=2, K=25){
  if (is.null(xmin)) xmin = min(x) - sd(x)
  if (is.null(xmax)) xmax = max(x) + sd(x)
  if (xmin > min(x)) xmin = min(x) - sd(x)
  if (xmax < max(x)) xmax = max(x) + sd(x)
  cnt = 0
  if (xmin < min(x)) cnt = cnt + 1
  if (xmax > max(x)) cnt = cnt + 1
  ##
  if (equid.knots){
    knots = seq(xmin,xmax,length=K-2)
    Dd = diff(diag(K),dif=pen.order) ## Penalty order
    Pd = t(Dd)%*%Dd #+ 1e-6*diag(KK)
  } else {
    knots = unique(c(xmin,quantile(x,p=seq(0,1,length=K-2-cnt)),xmax))
    ngrid = 501
    xgrid = seq(xmin,xmax,length=ngrid)
    d2B = D2Bsplines(xgrid,knots)
    Pd = 1/ngrid * (t(d2B)%*%d2B)
    pen.order = 2 ## Penalty based on 2nd derivatives of B-splines
  }
  ##
  ans = list(xmin=xmin,xmax=xmax,knots=knots,Pd=Pd,pen.order=pen.order)
  return(ans)
}
