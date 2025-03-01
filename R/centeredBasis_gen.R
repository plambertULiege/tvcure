#' Generation of a recentered cubic B-spline basis matrix.
#' @description Generation of a cubic B-spline basis matrix with recentered columns
#'  to handle the identifiability constraint in additive models. See Wood (CRC Press 2017, pp. 175-176) for more details.
#' @param x Vector of values where the "recentered" B-spline basis is evaluated.
#' @param knots Vector of knots that must cover the values in \code{x}.
#' @param cm (Optional) values subtracted from each column of the original B-spline matrix.
#' @param pen.order Order of the penalty applied on B-spline parameters. (Default: 2).
#'
#' @return List containing
#' \itemize{
#' \item \code{B} : centered cubic B-spline matrix obtained by subtracting \code{cm[j]} from the jth B-spline in column j of the original B-spline matrix evaluated at \code{x}.
#' \item \code{Dd} : difference matrix (of order \code{pen.order}) for the associated centered B-spline matrix.
#' \item \code{Pd} : penalty matrix (of order \code{pen.order}) for the associated centered B-spline matrix.
#' \item \code{K} : number of centered B-splines in the basis.
#' \item \code{cm} : values subtracted from each column of the original B-spline matrix. By default, this is a vector containing the mean of each column in the original B-spline matrix.
#'}
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#' @references Lambert, P. and Kreyenfeld, M. (2024). Exogenous time-varying covariates in double additive cure survival model
#' with application to fertility. \emph{Journal of the Royal Statistical Society, Series A}, under review.
#'
#' @examples
#' x = seq(0,1,by=.01)
#' knots = seq(0,1,length=5)
#' obj = centeredBasis.gen(x,knots)
#' matplot(x,obj$B,type="l",ylab="Centered B-splines")
#' colMeans(obj$B)
#'
#' @export
#'
centeredBasis.gen = function(x,knots,cm=NULL,pen.order=2){
  if ((max(x)>max(knots))|(min(x)<min(knots))){
    cat("The knots do no cover the values of x !!\n")
    return(NULL)
  }
  ##
  knots.x = knots ; pen.order.x = pen.order
  ##
  temp = Bsplines(x, knots.x)
  idx = 1:(ncol(temp)-1)
  B <- temp[,idx,drop=FALSE] ## Reduce the number of B-splines by 1 !!
  if (is.null(cm)){
    ## Consider a B-spline matrix over a dense equi-spaced set of values on range(knots)
    ##  to compute a recentering value per column of the matrix
    B.ref = Bsplines(seq(min(knots),max(knots),length=100), knots.x)[,idx]
    cm = colMeans(B.ref)
  }
  B = sweep(B,2,cm) ## Substract column mean from each column
  Dd = diff(diag(ncol(temp)),dif=pen.order.x)[,idx] ## Penalty order <---- Important to remove last column !!
  Pd = t(Dd) %*% Dd
  return(list(B=B,Dd=Dd,Pd=Pd,K=ncol(B),cm=cm))
}
