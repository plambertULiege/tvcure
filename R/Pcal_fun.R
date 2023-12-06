#' Function to generate a penalty matrix for additive terms.
#' @description Compute the penalty matrix associated to a vector containing fixed (non-penalized) parameters and equal-size sub-vectors of penalized spline parameters.
#' 
#' @param nfixed the number of fixed (i.e. non-penalized) parameters.
#' @param lambda a vector of \code{p} penalty parameters where each component is associated to a sub-vector of spline parameters of length \code{J}.
#' @param Pd.x a penalty matrix of size \code{J} associated to a given sub-vector of spline parameters.
#'.
#' @return A block diagonal penalty matrix of size \code{(nfixed+pJ)} given by Blockdiag(diag(0,\code{nfixed}), diag(\code{lambda}).kron.\code{Pd.x}).
#' 
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#' @references Lambert, P. and Kreyenfeld, M. (2024). Exogenous time-varying covariates in double additive cure survival model
#' with application to fertility. \emph{Journal of the Royal Statistical Society, Series A}, under review.
#' 
#' @examples
#' Dd = diff(diag(1,5),diff=2) ## Difference penalty matrix for a vector of length 5
#' Pd = t(Dd) %*% Dd ## Penalty matrix of order 2
#' nfixed = 2 ## 2 unpenalized parameters
#' ## Global penalty matrix when 2 unpenalized parameters and 2 additive terms with
#' ##   2 vectors of 5 P-splines coefficients with lambda values 10 and 100 respectively.
#' Pcal.fun(nfixed=2,lambda=c(10,100),Pd)
#' @export
#' 
Pcal.fun = function(nfixed,lambda,Pd.x){
  ntot = nfixed + ncol(Pd.x)*length(lambda)
  Pcal = matrix(0,ntot,ntot)
  Pcal[(nfixed+1):ntot,(nfixed+1):ntot] = kronecker(diag(lambda,length(lambda)),Pd.x)
  return(Pcal)
}
