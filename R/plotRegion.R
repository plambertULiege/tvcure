#' Generic function to plot a shaded region around values of a scalar function.
#' @usage plotRegion(x, mat,
#'            add=FALSE,xlim=range(x),ylim=range(mat),
#'            colfill="#D9D9D980",lwd=2,xlab="",ylab="",main="",...)

#' @param x n-vector with a grid of x values where the scalar function f(x) is evaluated.
#' @param mat (n x 3)-matrix containing on its ith row, the function value at x[i] and the bounds of an interval containing it, (f(x[i]), f.low(x[i]), f.up(x[i])).
#' @param add logical indicating if the shaded region should be superposed to an existing plot.
#' @param xlim x-limits. (Default: range of x).
#' @param ylim y-limits. (Default: range of mat).
#' @param colfill color used for filling the shaded region. (Default: "#D9D9D980").
#' @param lwd line width to plot (x,f(x)). (Default: 2).
#' @param xlab x-label. (Default: none).
#' @param ylab y-label. (Default: none).
#' @param main plot main title. (Default: none)
#' @param ... additional generic plotting arguments.
#'
#' @return No returned value (in addition to the plot).
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#' @references Lambert, P. and Kreyenfeld, M. (2025).
#' Time-varying exogenous covariates with frequently changing values in double additive cure survival model: an application to fertility.
#' \emph{Journal of the Royal Statistical Society, Series A}. <doi:10.1093/jrsssa/qnaf035>
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
#' obj = additive.tvcure(model) ## Extract additive terms
#'
#' ## Plot some of the fitted additive terms
#' ## par(mfrow=c(1,2))
#' with(obj$f1.grid$x1, plotRegion(x=x,mat=y.mat,xlab="x1",ylab="f(x1)"))
#' with(obj$f1.grid$x2, plotRegion(x=x,mat=y.mat,xlab="x2",ylab="f(x2)"))
plotRegion = function(x,mat,add=FALSE,xlim=range(x),ylim=range(mat),
                      colfill="#D9D9D980",lwd=2,xlab="",ylab="",main="",...){
  f = mat[,1] ; f.low = mat[,2] ; f.up = mat[,3]
  if (add==FALSE) plot(x,f,type="n",ylim=ylim,xlim=xlim,
                       lwd=lwd,xlab=xlab,ylab=ylab,main=main,...)
  polygon(c(x,rev(x)),c(f.low,rev(f.up)),col=colfill,border=F)
  ## adjustcolor("grey85", alpha.f=0.5): "#D9D9D980"
  ## adjustcolor("grey90", alpha.f=0.5): "#E5E5E580"
  lines(x,f,lwd=lwd)
}
