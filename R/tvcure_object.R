#' Object resulting from the fit of a tvcure model using function 'tvcure'.
#'
#' An object returned by the \code{\link{tvcure}} function: this is a list
#' with various components related to the fit of such a model.
#'
#' @return A \code{tvcure_object} is a list with following elements:
#' \itemize{
#' \item{\code{formula1} : \verb{ }}{A formula describing the linear predictor in the long-term (cure) survival (or quantum) submodel.}
#' \item{\code{formula2} : \verb{ }}{A formula describing the linear predictor in the short-term (cure) survival (or timing) submodel.}
#' \item{\code{method} : \verb{ }}{Method ("S0" or "F0") used to specify the dependence of the cumulative hazard dynamics on covariates.}
#' \item{\code{regr1} : \verb{ }}{List returned by \code{\link{DesignFormula}} when evaluated on \code{formula1}.}
#' \item{\code{regr2} : \verb{ }}{List returned by \code{\link{DesignFormula}} when evaluated on \code{formula2}.}
#' \item{\code{K0} : \verb{ }}{Number of B-splines used to specify \eqn{\log f_0(t)}.}
#' \item{\code{fit} : \verb{ }}{A list containing different elements describing the fitted tvcure model:}
#'   \itemize{
#'   \item{\code{llik} : \verb{ }}{Log likelihood value of the fitted tvcure model at convergence.}
#'   \item{\code{lpen} : \verb{ }}{Log of the penalized joint posterior at convergence.}
#'   \item{\code{dev} : \verb{ }}{Deviance of the fitted tvcure model at convergence.}
#'   \item{\code{res} : \verb{ }}{Standardized residual \eqn{(d_{ij}-\mu_{ij})/\sqrt{\mu_{ij}}} for unit \eqn{i} at time \eqn{t_{ij}} where \eqn{\mu_{j}=h_p(t_{ij}|z(t_{ij}),x(t_{ij}))} and \eqn{d_{ij}} is the event indicator.}
#'   \item{\code{phi} : \verb{ }}{Vector of length \eqn{K_0} containing the estimated B-splines coefficients in \eqn{\log f_0(t)}.}
#'   \item{\code{nbeta} : \verb{ }}{Number of regression and spline parameters in the long-term (cure) survival (or quantum) submodel.}
#'   \item{\code{beta} : \verb{ }}{(nbeta x 6) matrix containing the point estimates, standard errors, credible intervals, Z-scores and P-values of the regression and spline parameters in the long-term (cure) survival (or quantum) submodel.}
#'   \item{\code{ngamma} : \verb{ }}{Number of regression and spline parameters in the short-term (cure) survival (or timing) submodel.}
#'   \item{\code{gamma} : \verb{ }}{(ngamma x 6) matrix containing the point estimates, standard errors, credible intervals, Z-scores and P-values of the regression and spline parameters in the short-term (cure) survival (or timing) submodel.}
#'   \item{\code{gam} : \verb{ }}{ngamma-vector with the point estimates of the regression and spline parameters in the short-term (cure) survival (or timing) submodel.}
#'   \item{\code{grad.beta} : \verb{ }}{Gradient of the log joint posterior of <beta>, the regression and spline parameters in the long-term (cure) survival (or quantum) submodel.}
#'   \item{\code{Hes.beta} : \verb{ }}{Hessian of the log joint posterior of <beta>.}
#'   \item{\code{Hes.beta0} : \verb{ }}{Hessian of the log joint posterior of <beta> (with the roughness penalty part omitted).}
#' \item{\code{grad.gamma} : \verb{ }}{Gradient of the log joint posterior of <gamma>, the regression and spline parameters in the short-term (cure) survival (or timing) submodel.}
#'   \item{\code{Hes.gamma} : \verb{ }}{Hessian of the log joint posterior of <gamma>.}
#'   \item{\code{Hes.gamma0} : \verb{ }}{Hessian of the log joint posterior of <gamma> (with the roughness penalty part omitted).}
#'   \item{\code{Mcal.1} : \verb{ }}{Hessian of the log joint posterior of the spline parameters in <beta> conditionally on the non-penalized parameters.}
#'   \item{\code{Mcal.2} : \verb{ }}{Hessian of the log joint posterior of the spline parameters in <gamma> conditionally on the non-penalized parameters.}
#'   \item{\code{Hes.betgam} : \verb{ }}{(nbeta x ngamma) matrix with the cross derivatives of the log joint posterior of (<beta>,<gamma>).}
#'   \item{\code{grad.regr} : \verb{ }}{Gradient of the log joint posterior of <beta,gamma>.}
#'   \item{\code{Hes.regr} : \verb{ }}{Hessian of the log joint posterior of <beta,gamma>.}
#'   \item{\code{Hes.regr0} : \verb{ }}{Hessian of the log joint posterior of <beta,gamma> (with the roughness penalty part omitted).}
#'   \item{\code{grad.phi} : \verb{ }}{Gradient of the log joint posterior of <phi>, the spline parameters in \eqn{\log f_0(t)}.}
#'   \item{\code{Hes.phi} : \verb{ }}{Hessian of the log joint posterior of <phi>.}
#'   \item{\code{Hes.phi0} : \verb{ }}{Hessian of the log joint posterior of <phi> (with the roughness penalty part omitted).}
#'   \item{\code{T} : \verb{ }}{Follow-up time after which a unit is declared cured in the absence of a past event.}
#'   \item{\code{t.grid} : \verb{ }}{Grid of discrete time values on (1,T): 1,...,T.}
#'   \item{\code{f0.grid} : \verb{ }}{Estimated values for \eqn{f_0(t)} on \code{t.grid}.}
#'   \item{\code{F0.grid} : \verb{ }}{Estimated values for \eqn{F_0(t)} on \code{t.grid}.}
#'   \item{\code{S0.grid} : \verb{ }}{Estimated values for \eqn{S_0(t)} on \code{t.grid}.}
#'   \item{\code{dlf0.grid} : \verb{ }}{(ngrid x length(phi)) matrix with the jth line containing the gradient of \eqn{\log f_0(t_j)} w.r.t. <phi>.}
#'   \item{\code{dlF0.grid} : \verb{ }}{(ngrid x length(phi)) matrix with the jth line containing the gradient of \eqn{\log F_0(t_j)} w.r.t. <phi>.}
#'   \item{\code{dlS0.grid} : \verb{ }}{(ngrid x length(phi)) matrix with the jth line containing the gradient of \eqn{\log S_0(t_j)} w.r.t. <phi>.}
#'   \item{\code{k.ref} : \verb{ }}{Index of the reference component in <phi> set to 0.0.}
#'   \item{\code{a, b} : \verb{ }}{Hyperparameters of the Gamma(a,b) prior for the penalty parameters of the additive terms.}
#'   \item{\code{criterion} : \verb{ }}{Criterion used to assess convergence of the estimation procedure.}
#'   \item{\code{grad.psi} : \verb{ }}{Gradient of the log joint posterior of <phi[-k.ref]>, i.e. the spline parameters in \eqn{\log f_0(t)} with the fixed reference component omitted.}
#'   \item{\code{Hes.psi0} : \verb{ }}{Hessian of the log joint posterior of <phi[-k.ref]> (with the roughness penalty part omitted).}
#'   \item{\code{Hes.psi} : \verb{ }}{Hessian of the log joint posterior of <phi[-k.ref]>.}
#'   \item{\code{tau} : \verb{ }}{Selected value for the penalty parameter \eqn{\tau} tuning the smoothness of \eqn{\log f_0(t)}.}
#'   \item{\code{pen.order0} : \verb{ }}{Penalty order for the P-splines used to specify \eqn{\log f_0(t)}.}
#'   \item{\code{lambda1} : \verb{ }}{Selected values for the penalty parameters \eqn{\lambda_1} tuning the smoothness of the additive terms in the long-term (cure) survival (or quantum) submodel.}
#'   \item{\code{pen.order1} : \verb{ }}{Penalty order for the P-splines in the long-term survival (or quantum) submodel.}
#'   \item{\code{lambda2} : \verb{ }}{Selected values for the penalty parameters \eqn{\lambda_2} tuning the smoothness of the additive terms in the short-term (cure) survival (or timing) submodel.}
#'   \item{\code{pen.order2} : \verb{ }}{Penalty order for the P-splines in the short-term survival (or timing) submodel.}
#'   \item{\code{tau.method} : \verb{ }}{Method used to calculate the posterior mode of \eqn{p(\tau_0|{\cal D})}.}
#'   \item{\code{lambda.method} : \verb{ }}{Method used to select the penalty parameters of the additive terms in the long-term survival (or quantum) submodel.}
#'   \item{\code{ED1} : \verb{ }}{Effective degrees of freedom for each of the additive terms in the long-term survival (or quantum) submodel, with the selected statistical test for significance and its P-value.}
#'   \item{\code{ED2} : \verb{ }}{Effective degrees of freedom for each of the additive terms in the short-term survival (or timing) submodel, with the selected statistical test for significance and its P-value.}
#'   \item{\code{ED1.Tr} : \verb{ }}{Effective degrees of freedom for each of the additive terms in the long-term survival (or quantum) submodel, with Wood's statistical test for significance and its P-value.}
#'   \item{\code{ED2.Tr} : \verb{ }}{Effective degrees of freedom for each of the additive terms in the short-term survival (or timing) submodel, with Wood's statistical test for significance and its P-value.}
#'   \item{\code{ED1.Chi2} : \verb{ }}{Effective degrees of freedom for each of the additive terms in the long-term survival (or quantum) submodel, with a Chi-square test for significance and its P-value.}
#'   \item{\code{ED2.Chi2} : \verb{ }}{Effective degrees of freedom for each of the additive terms in the short-term survival (or timing) submodel, with a Chi-square test for significance and its P-value.}
#'   \item{\code{nobs} : \verb{ }}{Total number of observations.}
#'   \item{\code{n} : \verb{ }}{Total number of units or subjects.}
#'   \item{\code{d} : \verb{ }}{Total number of observed events.}
#'   \item{\code{ED1.tot} : \verb{ }}{Total effective degrees of freedom for the long-term survival (or quantum) submodel.}
#'   \item{\code{ED2.tot} : \verb{ }}{Total effective degrees of freedom for the short-term survival (or timing) submodel.}
#'   \item{\code{ED.tot} : \verb{ }}{Total effective degrees of freedom for the tvcure model.}
#'   \item{\code{AIC} : \verb{ }}{Akaike information criterion for the fitted model with a penalty calculated using the total effective degrees of freedom, -2log(L) + 2*ED.tot, larger values being preferred during model selection.}
#'   \item{\code{BIC} : \verb{ }}{Bayesian (Schwarz) information criterion for the fitted model with a penalty calculated using the total effective degrees of freedom and the total number of observed events, -2log(L) + log(d)*ED.tot, smaller values being preferred during model selection.}
#'   \item{\code{levidence} : \verb{ }}{Log-evidence of the fitted model, larger values being preferred during model selection.}
#'   \item{\code{iter} : \verb{ }}{Number of iterations required to achieve convergence.}
#'   \item{\code{elapsed.time} : \verb{ }}{Total duration (in seconds) of the estimation procedure.}
#'   }
#' \item{\code{call} : \verb{ }}{Function call.}
#' \item{\code{converged} : \verb{ }}{Binary convergence status.}
#' \item{\code{logLik} : \verb{ }}{Log-likelihood of the fitted model.}
#' }
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#' @references Lambert, P. and Kreyenfeld, M. (2024). Exogenous time-varying covariates in double additive cure survival model
#' with application to fertility. \emph{Journal of the Royal Statistical Society, Series A}, under review.
#'
#' @seealso \code{\link{tvcure}}, \code{\link{print.tvcure}}, \code{\link{plot.tvcure}}
#'
#' @name tvcure.object
NULL
