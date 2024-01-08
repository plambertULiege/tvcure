## require(cubicBsplines)
## require(numDeriv)
## require(ggplot2) ; require(gridExtra)
## require(Rfast)  ## Fast Crossprod
## require(Matrix) ## Sparse matrices
## require(tictoc) ## To evaluate computation times of code chunks:  tic("ici") ; ...code... ; toc()

## source("centeredBasis_gen.R")
## source("DesignFormula.R")
## source("s.R")
## source("qknots.R")
## source("Pcal_fun.R")
## source("additive.tvcure2.R")

## df.raw structure:
##  id time event     z1 z2 x1 x2
## NOTE: time must take positive integer values (1 being the unit of measurement for time)

#' Fit of a tvcure model.
#' @description Fit of a double additive cure survival model with exogenous time-varying covariates.
#' @usage tvcure(formula1, formula2, df,
#'        method=c("S0","F0"), K0=20, pen.order0=2,
#'        K1=10, pen.order1=2, K2=10, pen.order2=2,
#'        phi.0=NULL, beta.0=NULL, gamma.0=NULL,
#'        a.tau=1, b.tau=1e-6, a=1, b=1e-2,
#'        tau.0=100, tau.min=1, tau.method = c("LPS","LPS2","Schall","grid","none"),
#'        lambda1.0=NULL, lambda1.min=1, lambda2.0=NULL, lambda2.min=1,
#'        lambda.method=c("LPS","LPS2","LPS3","Schall","nlminb","none"),
#'        observed.hessian=TRUE, use.Rfast=TRUE, Wood.test=FALSE,
#'        ci.level=.95,
#'        criterion=c("levidence","deviance","lpen","AIC","BIC","gradient"),
#'        criterion.tol=1e-1, grad.tol=1e-2,
#'        iterlim=50, iter.verbose=TRUE, verbose=FALSE)
#' @param formula1 A formula describing the linear predictor in the long-term (cure) survival (or quantum) submodel.
#' @param formula2 A formula describing the linear predictor in the short-term (cure) survival (or timing) submodel.
#' @param df A data frame for survival data in a counting process format. It should always contain at least the following entries:
#' \itemize{
#' \item{\code{id} : \verb{ }}{the <id> of the unit associated to the data in a given line in the data frame.}
#' \item{\code{time} : \verb{ }}{the integer time at which the observations are reported. For a given unit, it should be a sequence of CONSECUTIVE integers starting at 1 for the first observation.}
#' \item{\code{event} : \verb{ }}{a sequence of 0-1 event indicators. For the lines corresponding to a given unit, it starts with 0 values concluded by a 0 in case of right-censoring or by a 1 if the event is observed at the end of the follow-up.}
#' }
#' @param method Method ("S0" or "F0") used to specify the dependence of the cumulative hazard dynamics on covariates (Default: "S0"):
#' Method S0: \eqn{S(t|x) = S_0(t)^{\exp^{\gamma'x}}} ;  Method F0: \eqn{F(t|x) = F_0(t)^{\exp{\gamma'x}}}
#' @param K0 Number of B-splines used to specify \eqn{\log f_0(t)} (Default: 20).
#' @param pen.order0 Penalty order for the P-splines used to specify \eqn{\log f_0(t)} (Default: 2).
#' @param K1 Number of P-splines for a given additive term in the long-term (or quantum) survival sumbodel (Default: 10).
#' @param pen.order1 Penalty order for the P-splines in the long-term survival (or quantum) sumbodel (Default: 2).
#' @param K2 Number of P-splines for a given additive term in the short-term (or timing) survival sumbodel (Default: 10).
#' @param pen.order2 Penalty order for the P-splines in the short-term survival (or timing) sumbodel (Default: 2).
#' @param phi.0 (Optional) vector of length \code{K0} with starting values for the P-spline parameters in \eqn{\log f_0(t)}.
#' @param beta.0 (Optional) starting value for the regression and spline parameters in the long-term survival (or quantum) submodel.
#' @param gamma.0 (Optional) starting value for the regression and spline parameters in the short-term survival (or timing) submodel.
#' @param a.tau Hyperprior parameter in the \eqn{Gamma(a.tau,b.tau)} prior for the penalty parameter \eqn{\tau} tuning the smoothness of \eqn{\log f_0(t)} (Default: 1.0).
#' @param b.tau Hyperprior parameter in the \eqn{Gamma(a.tau,b.tau)} prior for the penalty parameter \eqn{\tau} tuning the smoothness of \eqn{\log f_0(t)} (Default: 1e-6).
#' @param a Hyperprior parameter in the \eqn{Gamma(a,b)} priors for the penalty parameters \eqn{\lambda_1} and \eqn{\lambda_2} tuning the smoothness of the additive terms in the long-term (quantum) and short-term (timing) survival submodels. (Default: 1.0).
#' @param b Hyperprior parameter in the \eqn{Gamma(a,b)} priors for the penalty parameters \eqn{\lambda_1} and \eqn{\lambda_2} tuning the smoothness of the additive terms in the long-term (quantum) and short-term (timing) survival submodels. (Default: 1e-2).
#' @param tau.0 Starting value for \eqn{\tau}. (Default: 100).
#' @param tau.min Minimal value for the penalty parameter \eqn{\tau}. (Default: 1.0).
#' @param tau.method Method used to calculate the posterior mode of \eqn{p(\tau|data)}: "LPS", "LPS2", "Schall" (Fellner-Schall algorithm), "grid" (best choice in a regular grid on the log-scale) or "none" (stick to the initial value tau.0). LPS and LPS2, based on Laplace P-splines, both maximize the marginal posterior of the penalty parameter \eqn{\tau} using a fixed-point method, with LPS relying on the prior calculation of eigenvalues. (Default: "LPS").
#' @param lambda1.0 (Optional) J1-vector with starting values for the penalty parameters of the additive terms in the long-term survival (or quantum) submodel.
#' @param lambda1.min Minimal value for the J1 penalty parameters \eqn{\lambda_1} of the additive terms in the long-term survival (or quantum) submodel. (Default: 1.0).
#' @param lambda2.0 (Optional) J2-vector with starting values for the penalty parameters of the additive terms in the short-term survival (or timing) submodel.
#' @param lambda2.min Minimal value for the J2 penalty parameters \eqn{\lambda_2} of the additive terms in the short-term survival (or timing) submodel. (Default: 1.0).
#' @param lambda.method Method used ("LPS", "LPS2", "LPS3", "Schall", "nlminb" or "none") to select the penalty parameters of the additive terms in the long-term survival (or quantum) submodel:
#' \itemize{
#' \item{\code{LPS}, \code{LPS2}, or \code{LPS3} : \verb{ }}{based on Laplace P-splines where the marginal posterior of the penalty parameters is maximized using a fixed-point method. LPS is based on the prior calculation of eigenvalues (unlike LPS2) and delivers results of comparable quality to those of nlminb, but much more quickly. LPS3 work sequentially and separately on long- and short-term parameters with potentially convergence issues ;}
#' \item{\code{Schall} : \verb{ }}{Fellner-Schall method ;}
#' \item{\code{nlminb} : \verb{ }}{nonlinear maximization of the marginal posterior of the penalty parameters using the R function <nlminb> ;}
#' \item{\code{none} : \verb{ }}{penalty parameters are set at their initial values.}
#' }
#' @param observed.hessian Logical indicating if a fast approximation of the Hessian matrix based on cross-products is preferred over its expected value. (Default: TRUE).
#' @param use.Rfast Logical indicating if matrix functions from the Rfast package should be used to fasten computation. (Default: TRUE).
#' @param Wood.test Logical indicating if P-values based on Wood's test (Biometrika 2013) of the significance of additive terms should be preferred over basic Chi-square tests. (Default: FALSE).
#' @param ci.level Default value for the levels of the credible intervals. (Default: 0.95).
#' @param criterion Criterion used to assess convergence of the estimation procedure (Default: "levidence"):
#' \itemize{
#' \item{\code{levidence} : \verb{ }}{log of the evidence, i.e. of the marginal posterior of the penalty parameters at their selected values ;}
#' \item{\code{deviance} : \verb{ }}{deviance or -2 log(Likelihood) ;}
#' \item{\code{AIC} : \verb{ }}{Akaike information criterion ;}
#' \item{\code{BIC} : \verb{ }}{Bayesian (or Schwarz) information criterion ;}
#' \item{\code{gradient} : \verb{ }}{L2 norm of the gradient of the log of the joint posterior w.r.t. the regression and spline parameters.}
#' }
#' @param criterion.tol Maximum absolute difference between the successive values of the \code{criterion} values (when different from "gradient") to declare convergence. (Default: 1e-1).
#' @param grad.tol Tolerance value to declare convergence based on gradient values in an optimization procedure (such as Newton-Raphson). (Default: 1e-2).
#' @param iterlim Maximum number of iterations. (Default: 50).
#' @param iter.verbose Logical indicating if the values of the convergence criterions should be printed after each iteration. (Default: TRUE).
#' @param verbose Logical indicating if additional output based on gradients should be printed at the end of each iteration. (Default: FALSE).
#'
#' @return An object of type \code{\link{tvcure.object}}.
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#' @references Lambert, P. and Kreyenfeld, M. (2023). Exogenous time-varying covariates in double additive cure survival model
#' with application to fertility.
#' \emph{Journal of the Royal Statistical Society, Series A}, under review.
#'
#' @examples
#' require(tvcure)
#' ## Simulated data generation
#' beta = c(beta0=.4, beta1=-.2, beta2=.15) ; gam = c(gam1=.2, gam2=.2)
#' df.raw = simulateTVcureData(n=500, seed=123, beta=beta, gam=gam,
#'                           RC.dist="exponential",mu.cens=550)$df.raw
#' ## TVcure model fitting
#' tau.0 = 2.7 ; lambda1.0 = c(40,15) ; lambda2.0 = c(25,70) ## Optional
#' model = tvcure(~z1+z2+s(x1)+s(x2), ~z3+z4+s(x3)+s(x4), df=df.raw,
#'                tau.0=tau.0, lambda1.0=lambda1.0, lambda2.0=lambda2.0)
#' print(model)
#' plot(model, mfrow=c(1,2))
#'
#' @export
#'
tvcure = function(formula1, formula2, df,
                  method=c("S0","F0"),
                  K0=20, pen.order0=2,
                  K1=10, pen.order1=2,
                  K2=10, pen.order2=2,
                  phi.0=NULL, beta.0=NULL, gamma.0=NULL,
                  a.tau=1, b.tau=1e-6,  ## Prior on penalty parameter <tau> for log f0(t): Gamma(a.tau,b.tau)
                  a=1, b=1e-2, ## Prior on penalty parameters for the additive terms: Gamma(a,b)
                  tau.0=100, tau.min=1,
                  tau.method = c("LPS","LPS2","Schall","grid","none"),
                  lambda1.0=NULL, lambda1.min=1, lambda2.0=NULL, lambda2.min=1,
                  lambda.method=c("LPS","LPS2","LPS3","Schall","nlminb","none"), ## Penalty selection method for additive terms
                  observed.hessian=TRUE,  ## TRUE: X[event==1,]'X[event==1,] ; FALSE: X'diag(mu.ij)X
                  use.Rfast=TRUE,
                  Wood.test=FALSE,
                  ci.level=.95,
                  criterion=c("levidence","deviance","lpen","AIC","BIC","gradient"),
                  criterion.tol=1e-1,
                  grad.tol=1e-2,
                  iterlim=50,iter.verbose=TRUE,verbose=FALSE){
    ##
    cl <- match.call()
    aa = a ; bb = b
    method = match.arg(method) ## "S0": S(t|x) = S0(t)^exp(gamma'x) ;  "F0": F(t|x) = F0(t)^exp(gamma'x)
    criterion = match.arg(criterion) ## Criterion to assess convergence of the global algorithm
    tau.method = match.arg(tau.method)
    lambda.method = match.arg(lambda.method)
    if (missing(formula1)) {message("Missing model formula <formula1> for long-term survival !") ; return(NULL)}
    if (missing(formula2)) {message("Missing model formula <formula2> for short-term survival !") ; return(NULL)}
    ## if (attr(terms(formula2),"intercept")){
    ##     cat("Intercept automatically removed from <formula2> !\n")
    ##     formula2 = formula(paste(deparse(formula2),"-1")) ## Remove intercept from <formula2>
    ## }
    if (missing(df)) {message("Missing data frame <df> with <id>, <time>, event indicator <event>, and covariate values !") ; return(NULL)}
    ##
    if (!is.null(phi.0)) K0 = length(phi.0)
    if (!is.null(beta.0))  K1 = length(beta.0)
    if (!is.null(gamma.0)) K2 = length(gamma.0)
    ##
    id = unique(df$id) ## Subject ids
    n = length(id)     ## Number of subjects
    ## Make sure that <time> takes integer values !!
    temp = round(df$time)
    if (!all(temp == df$time)){
        stop("<time> must take integer values !!")
    }
    df$time = temp
    T = max(df$time) ## Maximum follow-up time (that should be RC)
    ## Significance level
    alpha = 1-ci.level
    z.alpha = qnorm(1-.5*alpha)
    ## Preliminary NP estimation of S(t)
    tab = with(df, table(time,event))
    n.risk = rowSums(tab) ; n.event = tab[,2]
    h0.hat = n.event / n.risk ; H0.hat = cumsum(h0.hat) ; lS0.hat = -H0.hat
    ## Regression models for long-term (formula1) & short-term (formula2) survival
    ## ---------------------------------------------------------------------------
    regr1 = DesignFormula(formula1, data=df, K=K1, pen.order=pen.order1) ##, n=n)
    regr2 = DesignFormula(formula2, data=df, K=K2, pen.order=pen.order2, nointercept=TRUE) ##, n=n)
    ## Extract key stuff from regr1 & regr2
    nobs = nrow(regr1$Xcal)     ## Number of data entries
    J1 = regr1$J ; J2 = regr2$J ## Number of additive terms
    K1 = regr1$K ; K2 = regr2$K ## Number of B-splines per additive term
    nfixed1 = regr1$nfixed ; nfixed2 = regr2$nfixed ## Number of 'non-spline' regression parameters
    lambda1.lab = regr1$lambda.lab ; lambda2.lab = regr2$lambda.lab ## Penalty labels
    Pd1.x = regr1$Pd.x ; Pd2.x = regr2$Pd.x ## Basic penalty matrix for an additive term
    Z1 = regr1$Z ; Z2 = regr2$Z ## Design matrices associated to 'non-spline' regression parameters
    regr1.lab = colnames(regr1$Xcal) ; regr2.lab = colnames(regr2$Xcal) ## Labels of the regression parms
    addregr1.lab = regr1$additive.lab ; addregr2.lab = regr2$additive.lab ## Labels of the additive terms
    ##
    ## Log-determinant of a positive definite matrix based on Choleski
    ##  Note: faster than methods based on functions svd, qr or determinant
    ## ---------------------------------------------------------------------
    ldet.fun = function(x) c(determinant(x,logarithm=TRUE)$mod)
    ## ldet.fun = function(x) 2*sum(log(diag(chol(x))))
    ## End ldet.fun
    ##
    ## Fast calculation of log |B'WB + lambda*Pd| using pre-computed eigenvalues:
    ##  log |B'WB + lambda*Pd| = ldet0 + sum(log(lambda1)) + sum_j(log(tau+dj))
    ## --------------------------------------------------------------------------

    ## Starting from B'WB and Pd
    ev.fun = function(BwB, Pd){  ## t(B) %*% diag(w) %*% B
        K = ncol(Pd) ; rk = qr(Pd)$rank ; id1 = 1:rk ; id2 = (1:K)[-id1]
        temp2 = svd(Pd) ; U = temp2$u ; lambda1 = temp2$d[id1]
        BwB = t(U) %*% (BwB %*% U)
        M = BwB[id1,id1] - BwB[id1,id2,drop=FALSE]%*%solve(BwB[id2,id2,drop=FALSE])%*%BwB[id2,id1,drop=FALSE]
        MM = sqrt(1/lambda1) * t(sqrt(1/lambda1)*t(M))
        ## MM2 = diag(1/sqrt(lambda1))%*%M%*%diag(1/sqrt(lambda1))
        dj = svd(MM)$d
        ldet0 = ldet.fun(BwB[id2,id2,drop=FALSE])
        ## log |B'WB + lambda*Pd| = ldet0 + sum(log(lambda1)) + sum_j(log(tau+dj))
        return(list(ldet0=ldet0,lambda1=lambda1,dj=dj))
    }
    ## Starting from B, W and Pd
    ev.fun2 = function(B, w=NULL, Pd){  ## t(B) %*% diag(w) %*% B
        if (is.null(w)) w = rep(1,nrow(B))
        K = ncol(Pd) ; rk = qr(Pd)$rank ; id1 = 1:rk ; id2 = (1:K)[-id1]
        temp2 = svd(Pd) ; U = temp2$u ; lambda1 = temp2$d[id1]
        Bt = (sqrt(w) * B) %*% temp2$u
        Bt1 = Bt[,id1] ; Bt0 = Bt[,-id1,drop=FALSE]
        B01 = t(Bt0)%*%Bt1
        M = t(Bt1)%*%Bt1 - t(B01)%*%solve(t(Bt0)%*%Bt0)%*%B01
        MM = sqrt(1/lambda1) * t(sqrt(1/lambda1)*t(M))
        ## MM2 = diag(1/sqrt(lambda1))%*%M%*%diag(1/sqrt(lambda1))
        dj = svd(MM)$d
        ldet0 = ldet.fun(t(Bt0)%*%Bt0)
        return(list(ldet0=ldet0,lambda1=lambda1,dj=dj))
    }
    ## Eigenvalues associated to additive terms
    id1 = as.logical(df$event) ##which(event==1)
    ev1.lst = ev2.lst = list()
    ## ... in long-term survival submodel
    if (J1 > 0){
        rk1 = qr(regr1$Pd.x)$rank
        ## Non-penalized Hessian for <beta>
        if (use.Rfast){
            Hes.beta0 = -Rfast::Crossprod(regr1$Xcal[id1,,drop=FALSE], regr1$Xcal[id1,,drop=FALSE])
        } else {
            Hes.beta0 = -crossprod(regr1$Xcal[id1,,drop=FALSE], regr1$Xcal[id1,,drop=FALSE])
        }
        for (j in 1:J1){ ## Loop over additive terms in the long-term survival sub-model
            idx = nfixed1 + (j-1)*K1 + (1:K1)
            Hes.betaj = Hes.beta0[idx,idx,drop=FALSE]
            ## if (use.Rfast){
            ##     Hes.betaj = -Rfast::Crossprod(regr1$Xcal[id1,idx,drop=FALSE], regr1$Xcal[id1,idx,drop=FALSE])
            ## } else {
            ##     Hes.betaj = -crossprod(regr1$Xcal[id1,idx,drop=FALSE], regr1$Xcal[id1,idx,drop=FALSE])
            ## }
            ev1.lst[[j]] = ev.fun(BwB=Hes.betaj,Pd=regr1$Pd.x)$dj
        }
    }
    ## ... in short-term survival submodel
    if (J2 > 0){
        rk2 = qr(regr2$Pd.x)$rank
        ## (Exact) Non-penalized Hessian for <gamma>:
        ## if (nogamma){ ## Absence of covariate (and no intercept for identification reasons)
        ##     eta.2 = 0.0
        ## } else {
        ##     eta.2 = c(regr2$Xcal %*% gamma) ## Linear predictors in short-term survival
        ## }
        ## if (method == "F0") tempo = exp(eta.2)*log(F0.grid[time])
        ## if (method == "S0") tempo = exp(eta.2)*log(S0.grid[time])
        ## XcalTempo = (1+tempo) * regr2$Xcal
        ## Hes.gamma0 = -Rfast::Crossprod(regr2$XcalTempo[id1,,drop=FALSE], regr2$XcalTemp[id1,,drop=FALSE])
        ##
        ## (Approximate) Non-penalized Hessian for <gamma>:
        if (use.Rfast){
            Hes.gamma0 = -Rfast::Crossprod(regr2$Xcal[id1,,drop=FALSE], regr2$Xcal[id1,,drop=FALSE])
        } else {
            Hes.gamma0 = -crossprod(regr2$Xcal[id1,,drop=FALSE], regr2$Xcal[id1,,drop=FALSE])
        }
        for (j in 1:J2){ ## Loop over additive terms in the long-term survival sub-model
            idx = nfixed2 + (j-1)*K2 + (1:K2)
            Hes.gamj = Hes.gamma0[idx,idx,drop=FALSE]
            ## if (use.Rfast){
            ##     Hes.gamj = -Rfast::Crossprod(regr2$Xcal[id1,idx,drop=FALSE], regr2$Xcal[id1,idx,drop=FALSE])
            ## } else {
            ##     Hes.gamj = -crossprod(regr2$Xcal[id1,idx,drop=FALSE], regr2$Xcal[id1,idx,drop=FALSE])
            ## }
            ev2.lst[[j]] = ev.fun(BwB=Hes.gamj,Pd=regr2$Pd.x)$dj
        }
    }
    ## Initial values
    ## --------------
    ## ... for the regression parameters
    q1 = ncol(regr1$Xcal) ## Total number of regression and spline parameters in long-term survival
    q2 = ncol(regr2$Xcal) ## Total number of regression and spline parameters in short-term survival
    if (is.null(q2)){
        nogamma = TRUE ## Check whether covariates are absent in <formula2> --> gamma=0
    } else {
        nogamma = FALSE
    }
    if (is.null(beta.0)){
        beta.0 = rep(0,q1)
        beta.0[1] = log(tail(H0.hat,1))
    }
    names(beta.0) = regr1.lab
    if (!nogamma){
        if (is.null(gamma.0)) gamma.0 = rep(0,q2)
        names(gamma.0) = regr2.lab
    }
    ## ... for the penalty parameters in regression part
    if ((J1 > 0) & is.null(lambda1.0)) lambda1.0 = rep(100,J1)
    if ((J2 > 0) & is.null(lambda2.0)) lambda2.0 = rep(100,J2)
    ##
    ## Pre-evaluated IB-splines basis for F0(t)
    ## ----------------------------------------
    ## t.grid = 1:T
    ## obj.knots = qknots(t.grid, xmin=0, xmax=max(t.grid), equid.knots = TRUE, pen.order=pen.order0, K=K0)
    obj.knots = qknots(1:T, xmin=0, xmax=T+1, equid.knots = TRUE, pen.order=pen.order0, K=K0)
    knots = obj.knots$knots ## Knots
    Pd = obj.knots$Pd ## Penalty matrix
    ##
    ## IB0.grid = IBsplines(t0=0,1:T,knots=knots) ## TxS matrix [Ib0_s(t)]
    ## IB.T = c(IBsplines(t0=0,max(knots),knots=knots))
    B0.grid = Bsplines(1:T,knots=knots) ## TxS matrix [b0_s(t)]
    if (is.null(phi.0)) phi.0 = rep(0,ncol(B0.grid)) ## c(solve(t(B0.grid)%*%B0.grid + Pd)%*%t(B0.grid)%*%dnorm(1:T,.5*T,T/7))
    ## if (is.null(phi.0)) phi.0 = c(solve(t(B0.grid)%*%B0.grid + Pd) %*% t(B0.grid) %*% log(h0.hat+1e-3))
    names(phi.0) = paste("phi",1:length(phi.0),sep="")
    k.ref = ceiling(K0/2)
    phi.0 = phi.0 - phi.0[k.ref]
    colnames(B0.grid) = names(phi.0)
    tB0B0 = t(B0.grid[df$time,-k.ref]) %*% B0.grid[df$time,-k.ref]
    ##
    ## -----------------------------------------------------------------------------------------
    ff = function(phi, beta, gamma,
                  tau, lambda1, lambda2,
                  Dbeta=FALSE, Dgamma=FALSE, Dphi=FALSE, D2phi=FALSE,
                  Dlambda=FALSE, hessian=TRUE){
    ## -----------------------------------------------------------------------------------------
        Delta = 1 ## Unit of measurement for df$time
        eta.grid = c(B0.grid %*% phi)
        pi.grid = exp(eta.grid) / sum(exp(eta.grid))
        B.tilde = t(t(B0.grid) - c(t(B0.grid)%*%pi.grid)) ## TxS matrix
        ##
        f0.grid = pi.grid / Delta
        F0.grid = cumsum(pi.grid)
        S0.grid = pmax(1e-6, 1 - F0.grid)
        ##
        PdPhi = c(Pd%*%phi)
        quad.phi = sum(phi*PdPhi)
        ##
        eta.1 = c(regr1$Xcal %*% beta)  ## Linear predictors in long-term survival
        if (nogamma){ ## Absence of covariate (and no intercept for identification reasons)
            eta.2 = 0.0
        } else {
            eta.2 = c(regr2$Xcal %*% gamma) ## Linear predictors in short-term survival
        }
        time = df$time ; event = df$event ; id1 = as.logical(event) ##which(event==1)
        ##
        switch(method,
            "F0" = {
                lhp = (eta.1 + eta.2) + ((exp(eta.2)-1.0) * log(F0.grid[time])) + log(f0.grid[time])
            },
            "S0" = {
                lhp = (eta.1 + eta.2) + ((exp(eta.2)-1.0) * log(S0.grid[time])) + log(f0.grid[time])
            }
        )
        ## Penalty matrices
        ## ----------------
        lambda1.cur = lambda1 ; lambda2.cur = lambda2
        if (J1 > 0) P1.cur = Pcal.fun(nfixed1,lambda1.cur,Pd1.x) ## Penalty matrix for location
        else P1.cur = 0
        ##
        if (J2 > 0) P2.cur = Pcal.fun(nfixed2,lambda2.cur,Pd2.x) ## Penalty matrix for dispersion
        else P2.cur = 0
        ## Log-likelihood
        ## --------------
        mu.ij = exp(lhp)*Delta
        llik = -sum(mu.ij) + sum(log(mu.ij[id1]))
        llik.y = -sum(event) + 0 ## ... as only 0-1 values for the Poisson counts --> y*log(y)=0
        dev = -2*(llik - llik.y)
        ## Penalized llik
        ## --------------
        lpen = llik
        lpen = lpen + (aa-1)*log(tau) - bb*tau ## Prior on <tau>
        lpen = lpen + .5*(K0-pen.order0)*log(tau) - .5*tau*quad.phi  ## Prior on (phi|tau)
        if (J1 > 0){
            ## Prior on <lambda1>
            lpen = lpen + sum(dgamma(lambda1.cur,aa,bb,log=TRUE)) ## Prior on <lambda1>
            ## Prior on (beta|lambda1)
            ev.beta = svd(P1.cur)$d
            lpen = lpen + .5*sum(log(ev.beta[ev.beta>1e-6])) - .5*sum(beta*c(P1.cur%*%beta))
        }
        if (J2 > 0){
            ## Prior on <lambda2>
            lpen = lpen + sum(dgamma(lambda2.cur,aa,bb,log=TRUE))
            ## Prior on (gamma|lambda2)
            ev.gam = svd(P2.cur)$d
            lpen = lpen + .5*sum(log(ev.gam[ev.gam>1e-6])) - .5*sum(gamma*c(P2.cur%*%gamma))
        }
        ##
        ## Derivatives of lpen wrt <beta>
        ## ------------------------------
        ## Gradient wrt <beta>
        grad.beta = NULL
        Hes.beta = Mcal.1 = NULL
        ## Hes.beta = Hes.beta0 = Mcal.1 = NULL
        if (Dbeta){
            if (use.Rfast){
                grad.beta  = Rfast::colsums(regr1$Xcal*(event-mu.ij))
            } else {
                grad.beta  = colSums(regr1$Xcal*(event-mu.ij))
            }
            ## Additive terms
            if (J1 > 0) grad.beta = grad.beta - c(P1.cur%*%beta)
        }
        ## Hessian wrt <beta>
        if (Dbeta & hessian){
            if (!observed.hessian){ ## Note: already computed when observed.hessian=TRUE
                W1 = mu.ij
                if (use.Rfast){
                    Hes.beta0 = -Rfast::Crossprod(regr1$Xcal, W1*regr1$Xcal)
                } else {
                    Hes.beta0 = -crossprod(regr1$Xcal, W1*regr1$Xcal)
                ## Hes.beta0 = -t(regr1$Xcal) %*% (W1*regr1$Xcal)
                }
            }
            ## if (observed.hessian){
            ##   if (use.Rfast){
            ##     Hes.beta0 = -Rfast::Crossprod(regr1$Xcal[id1,,drop=FALSE], regr1$Xcal[id1,,drop=FALSE]) ## ICI ICI
            ##   } else {
            ##       Hes.beta0 = -crossprod(regr1$Xcal[id1,,drop=FALSE], regr1$Xcal[id1,,drop=FALSE]) ## ICI ICI
            ##   }
            ## } else {
            ##     W1 = mu.ij
            ##     if (use.Rfast){
            ##         Hes.beta0 = -Rfast::Crossprod(regr1$Xcal, W1*regr1$Xcal)
            ##     } else {
            ##         Hes.beta0 = -crossprod(regr1$Xcal, W1*regr1$Xcal)
            ##     ## Hes.beta0 = -t(regr1$Xcal) %*% (W1*regr1$Xcal)
            ##     }
            ## }
            ## toc()
            ##
            ## Additive terms
            Hes.beta = Hes.beta0 ## Non-penalized Hessian Hes.beta0 for <beta>
            if (J1 > 0){
                Hes.beta = Hes.beta - P1.cur ## Penalized Hessian for <beta>
                ## Useful quantites for <lambda1> update
                if (Dlambda){
                    if (observed.hessian){
                      if (use.Rfast){
                          tZWB.1 = Rfast::Crossprod(Z1[id1,,drop=FALSE], regr1$Bcal[id1,,drop=FALSE])
                          tZWZ.1 = crossprod(Z1[id1,,drop=FALSE], Z1[id1,,drop=FALSE])
                      } else {
                          tZWB.1 = crossprod(Z1[id1,,drop=FALSE], regr1$Bcal[id1,,drop=FALSE])
                          tZWZ.1 = crossprod(Z1[id1,], Z1[id1,,drop=FALSE])
                      }
                    } else {
                      if (use.Rfast){
                          tZWB.1 = Rfast::Crossprod(Z1, W1*regr1$Bcal) ; tZWZ.1 = crossprod(Z1, W1*Z1)
                      } else {
                          tZWB.1 = crossprod(Z1, W1*regr1$Bcal) ; tZWZ.1 = crossprod(Z1, W1*Z1)
                          ## tZWB.1 = t(Z1)%*%(W1*regr1$Bcal) ; tZWZ.1 = t(Z1)%*%(W1*Z1)
                      }
                    }
                    inv.tZWZ.1 = try(solve(tZWZ.1),TRUE)
                    ## if (!is.matrix(inv.tZWZ.1)) return(stopthis())
                    ##
                    ## Mcal.1 required to update <lambda1>
                    if (observed.hessian){
                      if (use.Rfast){
                        Mcal.1 = Rfast::Crossprod(regr1$Bcal[id1,,drop=FALSE], regr1$Bcal[id1,,drop=FALSE])
                      } else {
                        Mcal.1 = crossprod(regr1$Bcal[id1,,drop=FALSE], regr1$Bcal[id1,,drop=FALSE])
                      }
                    } else {
                      if (use.Rfast){
                        Mcal.1 = Rfast::Crossprod(regr1$Bcal, W1*regr1$Bcal)
                      } else {
                          Mcal.1 = crossprod(regr1$Bcal, W1*regr1$Bcal)
                      }
                    }
                    if (use.Rfast){
                      Mcal.1 = Mcal.1 - Rfast::Crossprod(tZWB.1, inv.tZWZ.1%*%tZWB.1)
                    } else {
                        Mcal.1 = Mcal.1 - crossprod(tZWB.1, inv.tZWZ.1%*%tZWB.1)
                    }
                }
            }
        }
        ## Derivatives of lpen wrt <gamma>
        ## -------------------------------
        ## Gradient wrt <gamma>
        grad.gamma = NULL
        Hes.gamma = Hes.gamma0 = Mcal.2 = NULL
        if ((Dgamma) & (!nogamma)){
            switch(method,
                "F0" = {
                    tempo = exp(eta.2)*log(F0.grid[time])
                },
                "S0" = {
                    tempo = exp(eta.2)*log(S0.grid[time]) ## <-- Pack_v2
                }
                )
            XcalTempo = (1+tempo) * regr2$Xcal
            ##
            if (use.Rfast){
                grad.gamma = Rfast::colsums((event-mu.ij)*XcalTempo)
            } else {
                grad.gamma = colSums((event-mu.ij)*XcalTempo)
            }
            ## Additive terms
            if (J2 > 0) grad.gamma = grad.gamma - c(P2.cur%*%gamma) ## Penalized gradient for gamma
        }
        ## Hessian wrt <gamma>
        if ((Dgamma) & (!nogamma) & (hessian)){
            W2 = (event-mu.ij)*tempo - mu.ij*(1+tempo)^2
            if (observed.hessian){
              if (use.Rfast){
                Hes.gamma0 = -Rfast::Crossprod(XcalTempo[id1,,drop=FALSE], XcalTempo[id1,,drop=FALSE])
              } else {
                Hes.gamma0 = -crossprod(XcalTempo[id1,,drop=FALSE], XcalTempo[id1,,drop=FALSE])
              }
            } else {
              if (use.Rfast){
                Hes.gamma0 = Rfast::Crossprod(XcalTempo, W2*XcalTempo)
                ## Hes.gamma0 = Rfast::Crossprod(regr2$Xcal, W2*regr2$Xcal)
              } else {
                Hes.gamma0 = crossprod(XcalTempo, W2*XcalTempo)
                ## Hes.gamma0 = crossprod(regr2$Xcal, W2*regr2$Xcal)
              }
            }
            ## Additive terms
            Hes.gamma = Hes.gamma0 ## Non-penalized Hessian Hes.gamma0 for <gamma>
            if (J2 > 0){
                Hes.gamma = Hes.gamma - P2.cur ## Penalized Hessian for <gamma>
                ## Useful quantites for <lambda2> update
                if (Dlambda){
                    if (nfixed2 > 0){
                      if (use.Rfast){
                        tZWB.2 = Rfast::Crossprod(Z2, W2*regr2$Bcal) ; tZWZ.2 = Rfast::Crossprod(Z2, W2*Z2)
                      } else {
                        tZWB.2 = crossprod(Z2, W2*regr2$Bcal) ; tZWZ.2 = crossprod(Z2, W2*Z2)
                        ## tZWB.2 = t(Z2)%*%(W2*regr2$Bcal) ; tZWZ.2 = t(Z2)%*%(W2*Z2)
                      }
                      inv.tZWZ.2 = try(solve(tZWZ.2),TRUE)
                      ## if (!is.matrix(inv.tZWZ.2)) return(stopthis())
                      ##
                      if (use.Rfast){
                          Mcal.2 = Rfast::Crossprod(regr2$Bcal, W2*regr2$Bcal) - Rfast::Crossprod(tZWB.2, inv.tZWZ.2%*%tZWB.2)
                        } else {
                          Mcal.2 = crossprod(regr2$Bcal, W2*regr2$Bcal) - crossprod(tZWB.2  , inv.tZWZ.2%*%tZWB.2)
                          ## Mcal.2 = t(regr2$Bcal)%*%(W2*regr2$Bcal) -t(tZWB.2)%*%inv.tZWZ.2%*%tZWB.2
                        }
                    } else {
                      if (use.Rfast){
                        Mcal.2 = Rfast::Crossprod(regr2$Bcal, W2*regr2$Bcal)
                      } else {
                        Mcal.2 = crossprod(regr2$Bcal, W2*regr2$Bcal)
                      }
                    }
                }
            }
        }
        ## Cross-derivatives of lpen wrt <beta> & <gamma>
        ## ----------------------------------------------
        if (nogamma){
            grad.regr = grad.beta
            Hes.regr = Hes.beta
            Hes.regr0 = Hes.beta0
            Hes.betgam = NULL
        } else {
            grad.regr = c(grad.beta,grad.gamma)
            Hes.regr = Hes.regr0 = Hes.betgam = NULL
            if ((Dbeta) & (Dgamma) & (hessian)){
              if (use.Rfast){
                Hes.betgam = -Rfast::Crossprod(regr1$Xcal, (mu.ij*(1+tempo))*regr2$Xcal)
              } else {
                Hes.betgam = -crossprod(regr1$Xcal, (mu.ij*(1+tempo))*regr2$Xcal)
              }
              Hes.regr = rbind(cbind(Hes.beta, Hes.betgam),
                               cbind(t(Hes.betgam), Hes.gamma))
              Hes.regr0 = rbind(cbind(Hes.beta0, Hes.betgam),
                                cbind(t(Hes.betgam), Hes.gamma0))
              ## Hes.regr = rbind(cbind(Hes.beta, 0*Hes.betgam),
              ##                  cbind(0*t(Hes.betgam), Hes.gamma))
              ## Hes.regr0 = rbind(cbind(Hes.beta0, 0*Hes.betgam),
              ##                  cbind(0*t(Hes.betgam), Hes.gamma0))
            }
        }
        ## Derivatives of lpen wrt <phi>
        ## -----------------------------
        grad.phi = Hes.phi = Hes.phi0 = ed.phi = NULL
        dlf0.grid = dlF0.grid = dlS0.grid = NULL
        if (Dphi){
            S = K0
            ## Gradient
            dlf0.grid = B.tilde ## TxS matrix
            if (use.Rfast){
                temp = Rfast::colCumSums(pi.grid*B.tilde) ## TxS matrix
            } else {
                temp = apply(pi.grid*B.tilde,2,cumsum) ## TxS matrix
            }
            dlF0.grid =  temp / F0.grid ## TxS matrix
            dlS0.grid = -temp / S0.grid ## TxS matrix
            ##
            switch(method,
                   "F0" = {
                       dlhp = (exp(eta.2)-1.0) * dlF0.grid[time,] + dlf0.grid[time,] ## (nobs x S) matrix
                   },
                   "S0" = {
                       dlhp = (exp(eta.2)-1.0) * dlS0.grid[time,] + dlf0.grid[time,] ## (nobs x S) matrix
                   }
            )
            ## grad.phi = apply((event-mu.ij) * dlhp, 2,sum)
            ## grad.phi = colSums((event-mu.ij) * dlhp)
            ## toc() ; tic("  Dphi - P2")
            if (use.Rfast){
                grad.phi = Rfast::colsums(dlhp[id1,,drop=FALSE]) - Rfast::colsums(mu.ij * dlhp)
            } else {
                grad.phi = colSums(dlhp[id1,,drop=FALSE]) - colSums(mu.ij * dlhp)
            }
            grad.phi = grad.phi - tau * PdPhi ## Roughness penalty correction
            ## Hessian
            if (D2phi){
                ## Hessian without roughness penalty
                if (use.Rfast){
                    if (observed.hessian){
                        Hes.phi0 = - Rfast::Crossprod(dlhp[id1,,drop=FALSE], dlhp[id1,,drop=FALSE])
                    } else {
                        Hes.phi0 = - Rfast::Crossprod(dlhp, mu.ij * dlhp)
                    }
                } else {
                    if (observed.hessian){
                        Hes.phi0 = - crossprod(dlhp[id1,,drop=FALSE], dlhp[id1,,drop=FALSE])
                    } else {
                        Hes.phi0 = - crossprod(dlhp, mu.ij * dlhp)
                    }
                }
                ## Hessian with roughness penalty
                Hes.phi = Hes.phi0 - tau * Pd
                ## Effective dimension wrt <phi>
                ed.phi = sum(t(solve(-Hes.phi[-k.ref,-k.ref])) * (-Hes.phi0[-k.ref,-k.ref]))
            }
        }
        attr(phi,"ed.phi") = ed.phi ## Effective dimension for <phi> in F0(t) estimation
        ##
        ans = list(llik=llik, lpen=lpen, dev=dev, res=(event-mu.ij)/sqrt(mu.ij),
                   phi=phi, beta=beta, gamma=gamma,
                   nbeta=length(beta), ngamma=length(gamma),
                   grad.beta=grad.beta, Hes.beta=Hes.beta, Hes.beta0=Hes.beta0,
                   grad.gamma=grad.gamma, Hes.gamma=Hes.gamma, Hes.gamma0=Hes.gamma0,
                   Mcal.1=Mcal.1, Mcal.2=Mcal.2,
                   Hes.betgam=Hes.betgam,
                   grad.regr=grad.regr, Hes.regr=Hes.regr, Hes.regr0=Hes.regr0,
                   grad.phi=grad.phi, Hes.phi=Hes.phi, Hes.phi0=Hes.phi0,
                   T=T, t.grid=1:T, f0.grid=f0.grid, F0.grid=F0.grid, S0.grid=S0.grid,
                   dlf0.grid=dlf0.grid, dlF0.grid=dlF0.grid, dlS0.grid=dlS0.grid, k.ref=k.ref,
                   a=aa, b=bb)
        return(ans)
    } ## End ff
    ##
    ## Function 1: Penalty parameter selection using Laplace approximation to  p(lambda|data)
    ##    with an underlying computation of |B'WB + lamdba Pd| using precomputed eigenvalues
    ##    and the fixed point method to find lambda.MAP
    ## --------------------------------------------------------------------------------------
    select.lambda.LPS = function(itermax=50){
        ## Generic update of lambda
        update.lambda.fun = function(lambda,quad,ev,rk){
            ok.lam = FALSE ; iter.lam = 0
            while(!ok.lam){
                iter.lam = iter.lam + 1
                lambda.old = lambda
                ttr = sum(1 / (lambda+ev))
                lambda = (2*(aa-1) + rk) / (2*bb + quad + ttr)
                lam.dif = abs(lambda - lambda.old)
                ok.lam = (lam.dif < 1) | (iter.lam >= itermax)
            }
            return(lambda)
        }
        ## Update <lambda1> if (J1 > 0)
        if (J1 > 0){
            for (j in 1:J1){ ## Loop over additive terms in the long-term survival sub-model
                ## Update lambda1.cur
                idx = nfixed1 + (j-1)*K1 + (1:K1)
                theta.j = beta.cur[idx]
                quad.j = sum(theta.j*c(Pd1.x%*%theta.j))
                lambda1.cur[j] = update.lambda.fun(lambda1.cur[j],quad.j,ev1.lst[[j]],rk1)
            }
        }
        ## Update <lambda2> if (J2 > 0)
        if (J2 > 0){
            for (j in 1:J2){ ## Loop over additive terms in the short-term survival sub-model
                ## Update lambda2.cur
                idx = nfixed2 + (j-1)*K2 + (1:K2)
                theta.j = gamma.cur[idx]
                quad.j = sum(theta.j*c(Pd2.x%*%theta.j))
                lambda2.cur[j] = update.lambda.fun(lambda2.cur[j],quad.j,ev2.lst[[j]],rk2)
            }
        }
        ##
        return(list(lambda1=lambda1.cur,lambda2=lambda2.cur))
    } ## End select.lambda.LPS
    ##
    ## Function 2: Penalty parameter selection using Laplace approximation to  p(lambda|data)
    ##   with the fixed point method to find lambda.MAP (without precomputation of eigenvalues)
    ## ----------------------------------------------------------------------------------------
    select.lambda.LPS2 = function(Hes.regr0){
        ok.lam = FALSE ; iter.lam = 0
        ##
        update.Sigma = function(){
            if (J1 > 0) Pcal.1 = Pcal.fun(nfixed1,lambda1.cur,Pd1.x)
            else Pcal.1 = diag(0,nfixed1)
            if (J2 > 0) Pcal.2 = Pcal.fun(nfixed2,lambda2.cur,Pd2.x)
            else Pcal.2 = diag(0,nfixed2)
            ##
            ## Pcal = bdiag(Pcal.1,Pcal.2)
            ## Sigma = solve(-Hes.regr0 + Pcal)
            idx.1 = 1:ncol(Pcal.1)
            idx.2 = ncol(Pcal.1) + (1:ncol(Pcal.2))
            Sigma = as.matrix(bdiag(solve(-Hes.regr0[idx.1,idx.1] + Pcal.1), solve(-Hes.regr0[idx.2,idx.2] + Pcal.2)))
            ##
            return(Sigma)
        } ## End update.Sigma
        ##
        Sigma = update.Sigma() ## Update Sigma
        while(!ok.lam){
            iter.lam = iter.lam + 1
            ## Update <lambda1>
            lam1.dif = 0 ## Will be updated if (J1 > 0)
            if (J1 > 0){
                lambda1.old = lambda1.cur
                for (j in 1:J1){ ## Loop over additive terms in the long-term survival sub-model
                    ## Update lambda1.cur
                    idx = nfixed1 + (j-1)*K1 + (1:K1)
                    theta.j = beta.cur[idx]
                    quad.j = sum(theta.j*c(Pd1.x%*%theta.j))
                    ttr = max(sum(t(Sigma[idx,idx]) * Pd1.x), 1) ## max(sum(t(Sigma[idx,idx]) * Pd1.x),1)
                    ## Note: rank(Pd1.x) = (K1-pen.order1+1)
                    lambda1.cur[j] = ((aa-1) + K1-pen.order1+1) / (bb + quad.j + ttr)
                }
                lam1.dif = lambda1.cur-lambda1.old
            }
            ## Update <lambda2>
            lam2.dif = 0 ## Will be updated if (J2 > 0)
            if (J2 > 0){
                lambda2.old = lambda2.cur
                for (j in 1:J2){ ## Loop over additive terms in the short-term survival sub-model
                    ## Update lambda2.cur
                    idx = nfixed2 + (j-1)*K2 + (1:K2)
                    theta.j = gamma.cur[idx]
                    quad.j = sum(theta.j*c(Pd2.x%*%theta.j))
                    ttr = max(sum(t(Sigma[q1+idx, q1+idx]) * Pd2.x), 1) ## max(sum(t(Sigma[q1+idx, q1+idx]) * Pd2.x),1)
                    ## Note: rank(Pd2.x) = (K2-pen.order2+1)
                    lambda2.cur[j] = ((aa-1) + K2-pen.order2+1) / (bb + quad.j + ttr)
                }
                lam2.dif = lambda2.cur-lambda2.old
            }
            ## Update Sigma
            Sigma = update.Sigma()
            ##
            ## Stopping rule
            ok.lam = max(abs(c(lam1.dif,lam2.dif))) < 1 | (iter.lam > 20)
        } ## End while
        ##
        return(list(lambda1=lambda1.cur,lambda2=lambda2.cur,Sigma=Sigma,
                    niter=iter.lam,converged=ok.lam))
    } ## End select.lambda.LPS2
    ##
    ## Function 3: Penalty parameter selection using Laplace approximation to  p(lambda|data)
    ##     (separately for <lambda1> & <lambda2>
    ## --------------------------------------------------------------------------------------
    select.lambda.LPS3 = function(coef, nfixed, lambda, Pd, Mcal, pen.order, lambda.min, bb){
        Mcal.1 = Mcal
        lambda1.cur = lambda ; J1 = length(lambda1.cur)
        nfixed1 = nfixed
        Pd1.x = Pd ; K1 = ncol(Pd1.x)
        theta.cur = coef
        pen.order1 = pen.order
        lambda1.min = lambda.min
        ## Home-made Newton-Raphson
        ok.xi1 = FALSE ; iter.lam1 = 0
        while(!ok.xi1){
            iter.lam1 = iter.lam1 + 1
            P1.cur = Pcal.fun(nfixed1,lambda1.cur,Pd1.x) ## Update penalty matrix
            if (nfixed1 > 0){
                iMpen.1 = solve(Mcal.1 + P1.cur[-(1:nfixed1),-(1:nfixed1)])
            } else {
                iMpen.1 = solve(Mcal.1 + P1.cur)
            }
            ## Score.lambda1
            U.lam1 = U.xi1 = rep(0,J1)
            Rj.1 = list()
            ##
            for (j in 1:J1){
                lam1.j = rep(0,J1) ; lam1.j[j] = lambda1.cur[j]
                if (nfixed1 > 0){
                    Pj.1 = Pcal.fun(nfixed1,lam1.j,Pd1.x)[-(1:nfixed1),-(1:nfixed1)] ## Update penalty matrix
                } else {
                    Pj.1 = Pcal.fun(nfixed1,lam1.j,Pd1.x)
                }
                Rj.1[[j]] = iMpen.1%*%(Pj.1/lam1.j[j])
                idx = nfixed1 + (j-1)*K1 + (1:K1)
                theta.j = theta.cur[idx]
                quad1.cur = sum(theta.j*c(Pd1.x%*%theta.j)) + bb  ## <-----
                U.lam1[j] = .5*(K1-pen.order1+1)/lam1.j[j] -.5*quad1.cur -.5*sum(diag(Rj.1[[j]]))
            }
            ## Score.xi1  where  lambda1 = lambda1.min + exp(xi1)
            U.xi1 = U.lam1 * (lambda1.cur - lambda1.min)
            ##
            ## Hessian.lambda1
            Hes.lam1 = diag(0,J1)
            Only.diagHes1 = TRUE
            if (!Only.diagHes1){
                for (j in 1:J1){
                    for (k in j:J1){
                        temp = 0
                        if (j==k) temp = temp -.5*(K1-pen.order1+1)/lambda1.cur[j]^2
                        temp = temp + .5*sum(t(Rj.1[[j]])*Rj.1[[k]]) ## tr(XY) = sum(t(X)*Y)
                        Hes.lam1[j,k] = Hes.lam1[k,j] = temp
                    }
                }
            } else {
                for (j in 1:J1) Hes.lam1[j,j] = -.5*(K1-pen.order1+1)/lambda1.cur[j]^2 +.5*sum(t(Rj.1[[j]])*Rj.1[[j]])
            }
            ## Hessian.xi1  where  lambda1 = lambda1.min + exp(xi1)
            Hes.xi1 = t((lambda1.cur-lambda1.min) * t((lambda1.cur-lambda1.min) * Hes.lam1))
            diag(Hes.xi1) = diag(Hes.xi1) + U.lam1 * (lambda1.cur-lambda1.min)
            ## Update xi1  where  lambda1 = lambda1.min + exp(xi1)
            xi1.cur = log(lambda1.cur-lambda1.min) ## ; names(xi1.cur) = addregr1.lab ## paste("f.mu.",1:J1,sep="")
            dxi1 = -c(MASS::ginv(Hes.xi1) %*% U.xi1) ## Increment.xi1
            ##
            step = 1
            accept = FALSE
            while(!accept){
                xi1.prop = xi1.cur + step*dxi1
                accept = all(xi1.prop < 10) ## all(is.finite(exp(xi1.prop)))
                if (!accept) step = .5*step
            }
            if ((step !=1)&verbose) cat("Step-halving for lambda1: step=",step,"\n")
            xi1.cur = xi1.prop
            lambda1.cur = lambda1.min + exp(xi1.cur)
            ##
            ## Convergence ?
            ok.xi1 = (max(abs(U.xi1)) < grad.tol) | (iter.lam1 > 20)
        }
        ##
        return(list(lambda=lambda1.cur, U.lambda=U.lam1, Hes.lam=Hes.lam1,
                    xi=xi1.cur, U.xi=U.xi1, Hes.xi=Hes.xi1,
                    ok.xi=ok.xi1))
    } ## End select.lambda.LPS3
    ##
    ## Function testStat implements Wood (2013) Biometrika 100(1), 221-228
    ## Goal: evaluate H0: fj = Xt %*% beta = 0  when  (beta | D) ~ N(p,V)
    ## The type argument specifies the type of truncation to use.
    ## on entry `rank' should be an edf estimate
    ## 0. Default using the fractionally truncated pinv.
    ## 1. Round down to k if k<= rank < k+0.05, otherwise up.
    ## res.df is residual dof used to estimate scale. <=0 implies
    ## fixed scale.
    ## -------------------------------------------------------------------
    Tr.test = function(p,Xt,V,edf){
        ans <- testStat(p, Xt, V, min(ncol(Xt), edf), type = 0, res.df = -1)
        ## ans <- mgcv:::testStat(p, Xt, V, min(ncol(Xt), edf), type = 0, res.df = -1)
        return(ans)
    }
    ## End Tr.test
    ##
    ## Wood.test = function(p,Xt,V,edf){
    ##     ans <- testStat(p, Xt, V, min(ncol(Xt), edf), type = 0, res.df = -1)
    ##     return(ans)
    ## }
    ## End Wood.test
    ##
    ## Calculate EDF, Chi2 and Pval for additive terms
    ED.fun = function(fit, Wood.test=FALSE){
        ED1 = ED2 = ED1.Tr = ED2.Tr = ED1.Chi2 = ED2.Chi2 = NULL
        ##
        nbeta = fit$nbeta   ## Nbr of regression & spline parameters in long-term survival
        ngamma = fit$ngamma ## Nbr of regression & spline parameters in short-term survival
        ##
        Sigma.beta = with(fit, solve(-Hes.beta))
        ED.beta = rowSums(t(Sigma.beta) * (-fit$Hes.beta0))
        ##
        if (!nogamma){
            Sigma.gamma = with(fit, solve(-Hes.gamma))
            ED.gamma = rowSums(t(Sigma.gamma) * (-fit$Hes.gamma0))
        } else {
            Sigma.gamma = NULL
            ED.gamma = 0
        }
        ##
        ED.full = c(ED.beta,ED.gamma) ## Added on 2023.10.11
        ##
        ## Sigma.regr = with(fit, solve(-Hes.regr))
        ## ED.full = rowSums(t(Sigma.regr) * (-fit$Hes.regr0))
        ## Sigma.regr = with(fit, solve(-Hes.regr)) ## Removed on 2023.10.11
        ## ED.full = rowSums(t(Sigma.regr) * (-fit$Hes.regr0)) ## Removed on 2023.10.11
        ##
        if (J1 > 0){
            ## Sigma.beta = with(fit, solve(-Hes.beta))
            ## ED1.full = with(fit, rowSums(t(Sigma.beta) * (-Hes.beta0)))
            ## ED1.full = with(fit, diag(Sigma.beta %*% (-Hes.beta0)))
            ED1 = Chi2 = Tr = Pval.Tr = Pval.Chi2 = numeric(J1)
            for (j in 1:J1){
                idx = nfixed1 + (j-1)*K1 + (1:K1)
                beta.j = fit$beta[idx]
                ED1[j] = sum(ED.full[idx])
                ## Begin Wood
                ngrid = 200
                knots.x = regr1$knots.x[[j]] ; xL = min(knots.x) ; xU = max(knots.x)
                pen.order = regr1$pen.order
                x1.grid = seq(min(knots.x),max(knots.x),length=ngrid) ## Grid of values
                ## Centered B-spline basis
                cB = centeredBasis.gen(x1.grid,knots=knots.x,cm=NULL,pen.order)$B
                ## f1.grid[,j] = c(cB %*% beta.j)
                bele = Tr.test(beta.j,Xt=cB,Sigma.beta[idx,idx],ED1[j]) ## Added on 2023.10.11
                ## bele = Tr.test(beta.j,Xt=cB,Sigma.regr[idx,idx],ED1[j]) ## Removed on 2023.10.11
                Tr[j] = bele$stat
                Pval.Tr[j] = bele$pval
                ## cat("Wood:",Tr[j],Pval.Tr[j],ED1[j],"\n")
                ## End Wood
                ##
                ## "Naive" method
                Chi2[j] = sum(beta.j * c(solve(Sigma.beta[idx,idx]) %*% beta.j)) ## Added on 2023.10.11
                ## Chi2[j] = sum(beta.j * c(solve(Sigma.regr[idx,idx]) %*% beta.j)) ## Removed on 2023.10.11
                Pval.Chi2[j] = 1 - pchisq(Chi2[j],ED1[j])
                ## cat("Naive:",Chi2[j],Pval.Chi2[j],ED1[j],"\n")
                ## End Naive
            }
            ED1.Tr = cbind(ED1,Tr,Pval.Tr)
            ED1.Chi2 = cbind(ED1,Chi2,Pval.Chi2)
            rownames(ED1.Tr) = rownames(ED1.Chi2) = paste("f1(",regr1$additive.lab,")",sep="")
            colnames(ED1.Tr) = c("edf","Tr","Pval")
            colnames(ED1.Chi2) = c("edf","Chi2","Pval")
        }
        ##
        if (J2 > 0){
            ## Sigma.gamma = with(fit, solve(-Hes.gamma))
            ## ED2.full = with(fit, rowSums(t(Sigma.gamma) * (-Hes.gamma0)))
            ## ED2.full = with(fit, diag(Sigma.gamma %*% (-Hes.gamma0)))
            ED2 = Chi2 = Tr = Pval.Tr = Pval.Chi2 = numeric(J2)
            for (j in 1:J2){
                idx = nfixed2 + (j-1)*K2 + (1:K2)
                gamma.j = fit$gamma[idx]
                ED2[j] = sum(ED.full[nbeta + idx])
                ## Begin Wood
                ngrid = 200
                knots.x = regr2$knots.x[[j]] ; xL = min(knots.x) ; xU = max(knots.x)
                pen.order = regr2$pen.order
                x2.grid = seq(min(knots.x),max(knots.x),length=ngrid) ## Grid of values
                ## Centered B-spline basis
                cB = centeredBasis.gen(x2.grid,knots=knots.x,cm=NULL,pen.order)$B
                ## f2.grid[,j] = c(cB %*% gamma.j)
                bele = Tr.test(gamma.j,Xt=cB,Sigma.gamma[idx,idx],ED2[j]) ## Added on 2023.10.11
                ## bele = Tr.test(gamma.j,Xt=cB,Sigma.regr[nbeta+idx, nbeta+idx],ED2[j]) ## Removed on 2023.10.11
                Tr[j] = bele$stat
                Pval.Tr[j] = bele$pval
                ## cat("Wood:",Tr[j],Pval.Tr[j],ED2[j],"\n")
                ## End Wood
                ##
                ## "Naive" method
                Chi2[j] = sum(gamma.j * c(solve(Sigma.gamma[idx, idx]) %*% gamma.j)) ## Added on 2023.10.11
                ## Chi2[j] = sum(gamma.j * c(solve(Sigma.regr[nbeta+idx, nbeta+idx]) %*% gamma.j)) ## Removed on 2023.10.11
                Pval.Chi2[j] = 1 - pchisq(Chi2[j],ED2[j])
                ## cat("Naive:",Chi2[j],Pval[j],ED2[j],"\n")
                ## End Naive
            }
            ED2.Tr = cbind(ED2,Tr,Pval.Tr)
            ED2.Chi2 = cbind(ED2,Chi2,Pval.Chi2)
            rownames(ED2.Tr) = rownames(ED2.Chi2) = paste("f2(",regr2$additive.lab,")",sep="")
            colnames(ED2.Tr) = c("edf","Tr","Pval")
            colnames(ED2.Chi2) = c("edf","Chi2","Pval")
        }
        if (Wood.test){
            ED1=ED1.Tr ; ED2=ED2.Tr
        } else {
            ED1=ED1.Chi2 ; ED2=ED2.Chi2
        }
        ##
        return(list(ED1=ED1, ED2=ED2,
                    ED1.Tr=ED1.Tr, ED2.Tr=ED2.Tr,
                    ED1.Chi2=ED1.Chi2, ED2.Chi2=ED2.Chi2))
    } ## End ED.fun
    ##
    ## #############
    ## ESTIMATION ##
    ## #############
    L2norm = function(x) sqrt(sum(x^2))
    ## Generic Newton-Raphson algorithm
    ## --------------------------------
    NewtonRaphson = function(g, theta, tol=1e-2, itermax=15, verbose=FALSE){
        ntheta = length(theta)
        theta.cur = theta
        obj.cur = g(theta.cur,Dtheta=TRUE)
        g.start = obj.cur$g ## Function at the iteration start
        ## Convergence criterion using RDM (see e.g. Prague et al, 2013)
        RDM = with(obj.cur, sum(grad * dtheta)) / ntheta
        ok =  (RDM < tol^2) ## grad' (-H)^-1 grad / ntheta < tol^2 ?
        ## ok = (L2norm(obj.cur$grad) < tol) ## Stopping rule
        ## ok = all(abs(obj.cur$grad) < tol) ## Stopping rule
        iter = 0
        while(!ok){
            iter = iter + 1
            theta.cur = obj.cur$theta
            dtheta = c(obj.cur$dtheta)
            step = 1 ; nrep = 0
            repeat { ## Repeat step-halving directly till improve target function
                nrep = nrep + 1
                if (nrep > itermax) break ## if (nrep > 20) break
                theta.prop = theta.cur + step*dtheta ## Update.theta
                obj.prop = tryCatch(expr=g(theta.prop,Dtheta=TRUE), error=function(e) e)
                if (inherits(obj.prop, "error")){
                    step = .5*step
                } else {
                    if (obj.prop$g >= obj.cur$g) break
                    step = .5*step
                }
            }
            obj.cur = obj.prop
            ## Convergence criterion using RDM (see e.g. Prague et al, 2013)
            RDM = with(obj.cur, sum(grad * dtheta)) / ntheta
            ok =  (RDM < tol^2) ## grad' (-H)^-1 grad / ntheta < tol^2 ?
            ## Alternative stopping rules:
            ## ok = (L2norm(obj.cur$grad) < tol) ## Stopping rule
            ## ok = all(abs(obj.cur$grad) < tol) ## Stopping rule
            if (iter > itermax) break
        }
        if (verbose) cat(obj.cur$g," (niter = ",iter,") - grad = ",L2norm(obj.cur$grad),"\n",sep="")
        ans = list(val=obj.cur$g, val.start=g.start, theta=obj.cur$theta, grad=obj.cur$grad, iter=iter)
        return(ans)
    } ## End NewtonRaphson
    ##
    ## Starting values
    ## ---------------
    phi.cur = phi.0 ; beta.cur = beta.0 ; gamma.cur = gamma.0
    tau.cur = tau.0 ; lambda1.cur = lambda1.0 ; lambda2.cur = lambda2.0
    Hes.xi1 = Hes.xi2 = Hes.lam1 = Hes.lam2 = NULL
    tau.iter = c()
    ## Functions to handle identification issue in polytomial logistic estimation of <f0>
    ## ----------------------------------------------------------------------------------
    psi2phi = function(psi){
        ans = append(psi,values=0,after=k.ref-1)
        names(ans) = names(phi.0)
        return(ans)
    }
    phi2psi = function(phi) return((phi-phi[k.ref])[-k.ref])
    ## ## Check whether there is no covariate in the short-term survival part
    ## ##  in which case  <gamma> should be a scalar fixed to 0
    ## nogamma = (regr2$nfixed==1) & (colnames(regr2$Xcal)[1]=="(Intercept)") & (regr2$J == 0)
    ##
    converged = FALSE ; iter = 0
    final.iteration = FALSE ## Idea: when the EDs of the additive terms stabilize,
    ##
    ptm <- proc.time() ## Start timer
    ##
    if (tau.method == "grid") tau.stable = FALSE ## Initiate the algorithm with the selection of <tau> for F0(t)
    tau.iter = tau.cur
    dev.old = AIC.old = BIC.old = 1e12 ## Necessary for convergence criterion based on <deviance>
    lpen.old = levidence.old = -1e12    ##  ... or <AIC> or <BIC> or <lpen> or <log(evidence)>
    ## cat("Start of the iterative procedure\n")
    while(!converged){ ## Global estimation loop
        iter = iter + 1
        ## =====================
        ## -1- Estimation of F0  (through <phi> and <psi>)
        ## =====================
        ## cat("phi\n")
        ## tic("1: phi:")
        nlm.version = FALSE
        ##
        if (nlm.version){
            ## <nlm> estimation of <phi> (and <psi>) (F0 estimation)
            ## -------------------------------------
            g.psi.nlm = function(psi, tau){
                phi = psi2phi(psi)
                obj.cur = ff(phi=phi, beta=beta.cur, gamma=gamma.cur,
                             tau=tau,
                             lambda1=lambda1.cur, lambda2=lambda2.cur,
                             Dphi=TRUE,D2phi=FALSE)
                ans = -obj.cur$lpen
                attr(ans, "gradient") = -obj.cur$grad.phi[-k.ref]
                return(ans)
            } ## End g.psi.nlm
            psi.cur = phi2psi(phi.cur)
            psi.nlm = nlm(f=g.psi.nlm,psi.cur,tau=tau.cur) ##,hessian=TRUE) ## <---------
            psi.cur = psi.nlm$est ; phi.cur = psi2phi(psi.cur)
        }
        ##
        if (!nlm.version){
            ## Newton-Raphson for <phi> (F0 estimation)
            ## ------------------------
            g.psi = function(theta,Dtheta=TRUE){
                phi = psi2phi(theta)
                obj.cur = ff(phi=phi, beta=beta.cur, gamma=gamma.cur,
                             tau=tau.cur, lambda1=lambda1.cur, lambda2=lambda2.cur,
                             Dphi=Dtheta, D2phi=Dtheta, Dbeta=FALSE, Dgamma=FALSE)
                grad = Sigma = dtheta = NULL
                if (Dtheta){
                    grad = obj.cur$grad.phi[-k.ref]
                    ## A = Matrix::nearPD(-obj.cur$Hes.phi[-k.ref,-k.ref])$mat
                    ## dtheta = solve(A, grad)
                    A = -obj.cur$Hes.phi[-k.ref,-k.ref] + diag(1e-6,length(grad))
                    dtheta = solve(A, grad)
                    attr(theta,"ed.phi") = attr(obj.cur$phi,"ed.phi")
                }
                ans = list(g=obj.cur$lpen, theta=theta, dtheta=dtheta, grad=grad)
                return(ans)
            } ## End g.psi
            ##
            psi.cur = phi2psi(phi.cur)
            psi.NR = NewtonRaphson(g=g.psi,theta=psi.cur,tol=grad.tol)
            psi.cur = psi.NR$theta ; phi.cur = psi2phi(psi.cur)
            ## ed.phi = attr(psi.cur,"ed.phi") ## Effective dim of <phi> in the estimation of F0(t)
        }
        ##
        ## Selection of <tau> (Penalty parameter to estimate <phi> in F0)
        ## ------------------
        switch(tau.method,
               "LPS" = {
                   ## Generic update of <tau>
                   update.tau.fun = function(tau,quad,ev,rk,itermax=50){
                       ok = FALSE ; iter = 0
                       while(!ok){
                           iter = iter + 1
                           tau.old = tau
                           ttr = sum(1 / (tau+ev))
                           tau = (2*(a.tau-1) + rk) / (2*b.tau + quad + ttr)
                           tau.dif = abs(tau - tau.old)
                           ok = (tau.dif < 1) | (iter >= itermax)
                       }
                       return(tau)
                   }
                   ## Update tau.cur
                   obj.cur = ff(phi=phi.cur, beta=beta.cur, gamma=gamma.cur,
                                tau=tau.cur, lambda1=lambda1.cur, lambda2=lambda2.cur,
                                Dphi=TRUE, D2phi=TRUE, Dbeta=FALSE, Dgamma=FALSE)
                   BwB = -obj.cur$Hes.phi0[-k.ref,-k.ref]
                   psi.cur =  phi2psi(phi.cur) ; quad = sum(psi.cur * c(Pd[-k.ref,-k.ref] %*% psi.cur))
                   rk0 = qr(Pd[-k.ref,-k.ref])$rank ## Rank of Penalty matrix
                   ev0 = ev.fun(BwB=BwB,Pd=Pd[-k.ref,-k.ref])$dj ## Eigenvalues for update of <tau>
                   tau.cur = update.tau.fun(tau.cur,quad,ev0,rk0) ## <tau> udpate
               },
               "LPS2" = {
                   obj.cur = ff(phi=phi.cur, beta=beta.cur, gamma=gamma.cur,
                                tau=tau.cur, lambda1=lambda1.cur, lambda2=lambda2.cur,
                                Dphi=TRUE, D2phi=TRUE, Dbeta=FALSE, Dgamma=FALSE)
                   BWB = -obj.cur$Hes.phi0[-k.ref,-k.ref]
                   ok.tau = FALSE
                   while(!ok.tau){
                       Sigma = solve(BWB+tau.cur*Pd[-k.ref,-k.ref])
                       ttr = sum(t(Sigma)*Pd[-k.ref,-k.ref])
                       quad = sum(phi.cur * c(Pd%*%phi.cur))
                       tau.old = tau.cur
                       tau.cur = max(tau.min, (K0-pen.order0) / (quad + ttr))
                       ## cat("tau:",tau.old," --> ",tau.cur,"\n")
                       ok.tau = (abs(tau.old-tau.cur) < 1)
                       ## if(ok.tau) cat("\n")
                   }
               },
               "Schall" = {
                   obj.cur = ff(phi=phi.cur, beta=beta.cur, gamma=gamma.cur,
                                tau=tau.cur, lambda1=lambda1.cur, lambda2=lambda2.cur,
                                Dphi=TRUE, D2phi=TRUE, Dbeta=FALSE, Dgamma=FALSE)
                   BWB = -obj.cur$Hes.phi0[-k.ref,-k.ref]
                   ed.phi = sum(t(solve(BWB+tau.cur*Pd[-k.ref,-k.ref]))*BWB) ## = sum(diag(solve(BWB+tau*BWB)%*%Pd))
                   quad = sum(phi.cur * c(Pd%*%phi.cur))
                   tau.cur = max(tau.min, (ed.phi-pen.order0)/quad)
               },
               "grid" = {
                   if (!tau.stable){
                       ## Starting grid for <tau> = roughness penalty parameter in F0(t)
                       tau.ref = tau.cur
                       tau.grid = 2^c((log2(tau.cur)-1):min(15,(log2(tau.cur)+1)))
                       ngrid = length(tau.grid)
                       ED.grid = BIC.grid = AIC.grid = numeric(ngrid)
                       psi.grid = c()
                       psi.cur = phi2psi(phi.cur)
                       for (j in 1:ngrid){
                           ## N-R
                           tau.cur = tau.grid[j]
                           if (tau.cur == tau.ref){
                               psi.NR2 = psi.NR
                           } else {
                               psi.NR2 = NewtonRaphson(g=g.psi,theta=psi.cur,tol=grad.tol)
                           }
                           psi.prop = psi.NR2$theta ; phi.prop = psi2phi(psi.prop)
                           psi.grid[[j]] = psi.prop
                           obj = ff(phi=phi.prop, beta=beta.cur, gamma=gamma.cur,
                                    tau=tau.grid[j], lambda1=lambda1.cur, lambda2=lambda2.cur,
                                    Dphi=TRUE,D2phi=TRUE)
                           ED.grid[j] = sum(diag(solve(obj$Hes.phi0[-k.ref,-k.ref] - tau.grid[j]*Pd[-k.ref,-k.ref]) %*% obj$Hes.phi0[-k.ref,-k.ref]))
                           BIC.grid[j] = obj$dev + ED.grid[j]*log(sum(n.event))
                           AIC.grid[j] = obj$dev + 2*ED.grid[j]
                       }
                       ##
                       jj = which.min(AIC.grid)
                       tau.cur = tau.grid[jj]
                       tau.iter = c(tau.iter, tau.cur)
                       tau.stable = (sum(diff(tail(tau.iter,4)) == 0) == 2) ## Stop selecting tau when unchanged 3 iterations in a row
                       psi.cur = psi.grid[[jj]] ; phi.cur = psi2phi(psi.cur)
                   }
               },
               "none"={
                   ## <tau.cur> unchanged
               }
        ) ## End of switch(tau.method)
        ##
        ## =======================================
        ## -2- Newton-Raphson for <beta> & <gamma>
        ## =======================================
        ## cat("beta - gamma\n")
        ## toc() ; tic("2: regr")
        idx1 = 1:length(beta.0)
        ##
        ## ## Select <beta,gamma> using nlminb and (gradient + Hessian) (function ff)
        ## ## -----------------------------------------------------------------------
        ## gg = function(theta){
        ##     beta = theta[idx1] ; gamma = theta[-idx1]
        ##     -ff(phi=phi.cur, beta=beta, gamma=gamma,tau=tau.cur, lambda1=lambda1.cur, lambda2=lambda2.cur,Dbeta=FALSE,Dgamma=FALSE)$lpen
        ## }
        ## gr = function(theta){
        ##     beta = theta[idx1] ; gamma = theta[-idx1]
        ##     -ff(phi=phi.cur, beta=beta, gamma=gamma,tau=tau.cur, lambda1=lambda1.cur, lambda2=lambda2.cur,Dbeta=TRUE,Dgamma=TRUE)$grad.regr
        ## }
        ## gh = function(theta){
        ##     beta = theta[idx1] ; gamma = theta[-idx1]
        ##     -ff(phi=phi.cur, beta=beta, gamma=gamma,tau=tau.cur, lambda1=lambda1.cur, lambda2=lambda2.cur,Dbeta=TRUE,Dgamma=TRUE)$Hes.regr
        ## }
        ## bele = nlminb(start=c(beta.cur,gamma.cur), objective=gg, gradient=gr, hessian=gh, control=list(eval.max=500,iter.max=500))
        ## beta.cur = bele$par[idx1] ; gamma.cur = bele$par[-idx1]
        ##
        ## Select <beta,gamma> using homemade Newton-Raphson
        ## --------------------------------------------------
        g.regr = function(theta,Dtheta=TRUE){
            beta = theta[idx1] ; gamma = theta[-idx1]
            obj.cur = ff(phi=phi.cur, beta=beta, gamma=gamma,
                      tau=tau.cur, lambda1=lambda1.cur, lambda2=lambda2.cur,
                      Dbeta=Dtheta,Dgamma=Dtheta)
            ## grad = obj.cur$grad.regr ; Sigma = diag(diag(-solve(obj.cur$Hes.regr)))
            ## grad = obj.cur$grad.regr ; Sigma = -solve(obj.cur$Hes.regr-diag(1e-6,length(theta)))
            grad = dtheta = NULL
            if (Dtheta){
                grad = obj.cur$grad.regr
                ## A = Matrix::nearPD(-obj.cur$Hes.regr)$mat
                A = -obj.cur$Hes.regr
                dtheta = solve(A, grad)
            }
            ##
            ans = list(g=obj.cur$lpen, theta=theta, dtheta=dtheta, grad=grad)
            return(ans)
        } ## End g.regr
        ## cat("Regr estimation: ")
        regr.NR = NewtonRaphson(g=g.regr,theta=c(beta.cur,gamma.cur),tol=grad.tol)
        ## cat("done in ",regr.NR$iter," iterations\n",sep="")
        beta.cur = regr.NR$theta[idx1] ; gamma.cur = regr.NR$theta[-idx1]
        ##
        ## ================================
        ## -3- Update <lambda1> & <lambda2> (i.e. penalty vectors for additive terms in short- and long-term survival)
        ## ================================
        ##
        ## cat("lambda\n")
        ## toc() ; tic("3: lambda")
        itermin = 1 ## Only start updating after <itermin> iterations
        update.lambda = ifelse(iter <= itermin, FALSE, (J1 > 0)|(J2 > 0))
        ##
        ## Evaluate necessary objects for the selection of the penalty parameters
        obj.cur = ff(phi.cur, beta.cur, gamma.cur,
                     tau=tau.cur, lambda1=lambda1.cur, lambda2=lambda2.cur,
                     Dphi=FALSE, Dbeta=TRUE, Dgamma=!nogamma,
                     Dlambda=update.lambda & ((lambda.method == "LPS3")))
        ##           Dlambda=update.lambda & ((lambda.method == "LPS3")|(lambda.method == "LPS2")))
        ##
        ## Method 1: LPS
        ## -------------
        if (lambda.method == "LPS"){ ## Laplace's method with fast determinant computation
            if (update.lambda & !final.iteration){
                ## -3- lambda1 & lambda2
                ## ---------------------
                if (J2 >0){
                    Hes.gamma0 = obj.cur$Hes.gamma0 ## Update Hes.gamma0 value
                    for (j in 1:J2){ ## Loop over additive terms in the long-term survival sub-model
                        idx = nfixed2 + (j-1)*K2 + (1:K2)
                        Hes.gamj = Hes.gamma0[idx,idx,drop=FALSE]
                        ## Recompute eigenvalues for update of lambda2 (those for lambda1 remain unchanged!)
                        ev2.lst[[j]] = ev.fun(BwB=Hes.gamj,Pd=regr2$Pd.x)$dj
                    }
                }
                if (J1 > 0 | J2 > 0){
                    temp = select.lambda.LPS()
                    lambda1.cur = temp$lambda1
                    lambda2.cur = temp$lambda2
                }
            }
        } ## Endif lambda.method == "LPS"
        ##
        ## Method 2: LPS2
        ## --------------
        if (lambda.method == "LPS2"){ ## Laplace's method (jointly for (lambd1,lambda2))
            if (update.lambda & !final.iteration){
                ## -3- lambda1 & lambda2
                ## ---------------------
                if (J1 > 0 | J2 > 0){
                    temp = select.lambda.LPS2(Hes.regr0=obj.cur$Hes.regr0)
                    lambda1.cur = temp$lambda1
                    lambda2.cur = temp$lambda2
                }
            }
        } ## Endif lambda.method == "LPS2"
        ##
        ## Method 3: LPS3
        ## ---------------
        if (lambda.method == "LPS3"){ ## Laplace's method (separately for <lambda1> & <lambda2>
            if (update.lambda & !final.iteration){
                ## -3a- lambda1 (long-term survival)
                ## ------------
                if (J1 > 0){
                    temp = select.lambda.LPS3(coef=beta.cur, nfixed=nfixed1,
                                             lambda=lambda1.cur, Pd=Pd1.x, Mcal=obj.cur$Mcal.1,
                                             pen.order=pen.order1, lambda.min=lambda1.min, bb=bb)
                    lambda1.cur = temp$lambda ## Vector of penalty parameters for the additive terms
                    Hes.xi1 = temp$Hes.xi
                    Hes.lam1 = temp$Hes.lam
                }
                ## -3b- lambda2 (short-term survival)
                ## ------------
                if (J2 > 0){
                   temp = select.lambda.LPS3(coef=gamma.cur, nfixed=nfixed2,
                                             lambda=lambda2.cur, Pd=Pd2.x, Mcal=obj.cur$Mcal.2,
                                             pen.order=pen.order2, lambda.min=lambda2.min, bb=bb)
                   lambda2.cur = temp$lambda ## Vector of penalty parameters for the additive terms
                   Hes.xi2 = temp$Hes.xi
                   Hes.lam2 = temp$Hes.lam
                }

            }
        } ## Endif lambda.method == "LPS3"
        ##
        ##
        ## Method 2: SCHALL
        ## ----------------
        if (lambda.method == "Schall"){ ## Schall's method
            ##
            ED.cur = ED.fun(obj.cur,Wood.test=Wood.test) ## Evaluate current ED for additive terms
            ED1.tot = regr1$nfixed + ifelse(J1==0, 0, sum(ED.cur$ED1[,1]))
            ED2.tot = ifelse(nogamma, 0, regr2$nfixed) + ifelse(J2==0, 0, sum(ED.cur$ED2[,1]))
            ED.tot = ED1.tot + ED2.tot
            ##
            if (update.lambda & !final.iteration){
                ## -3a- lambda1 (long-term survival)
                ## ------------
                if (J1 > 0){
                    for (j in 1:J1){
                        idx = nfixed1 + (j-1)*K1 + (1:K1)
                        beta.j = beta.cur[idx] ## Centered B-splines coefs for jth additive term
                        quad = sum(beta.j * c(Pd1.x %*% beta.j))
                        ed = ED.cur$ED1[j,1]
                        tau2 = quad / (ed-pen.order1)
                        sigma2 = sum(obj.cur$res^2) / (length(obj.cur$res)-ED.tot)
                        ## sigma2 = obj.cur$dev / (length(obj.cur$mu.ij)-ED.tot)
                        lambda1.cur[j] = max(lambda1.min, sigma2/tau2)
                    }
                }
                ## -3b- lambda2 (short-term survival)
                ## ------------
                if (J2 > 0){
                    for (j in 1:J2){
                        idx = nfixed2 + (j-1)*K2 + (1:K2)
                        gamma.j = gamma.cur[idx] ## Centered B-splines coefs for jth additive term
                        quad = sum(gamma.j* c(Pd2.x %*% gamma.j))
                        ed = ED.cur$ED2[j,1]
                        tau2 = quad / (ed-pen.order2)
                        sigma2 = sum(obj.cur$res^2) / (length(obj.cur$res)-ED.tot)
                        ## sigma2 = obj.cur$dev / (length(obj.cur$mu.ij)-ED.tot)
                        lambda2.cur[j] = max(lambda2.min, sigma2/tau2)
                    }
                }
            } ## Endif update.lambda
        } ## Endif lambda.method == "Schall"
        ##
        ## Method 3: nlminb
        ## ----------------
        if (lambda.method == "nlminb"){ ## Direct optimization of log-evidence
            ## Loss function to select <lambda>:
            ##     loss.fn(nu) = -log p(lambda=exp(nu) | data)
            loglambda.loss = function(loglambda,data){
                if (J1 > 0) lambda1 = exp(loglambda[1:J1])
                if (J2 > 0) lambda2 = exp(loglambda[(J1+1):(J1+J2)])
                obj.cur = ff(phi.cur, beta.cur, gamma.cur,
                             tau=tau.cur, lambda1=lambda1, lambda2=lambda2,
                             Dphi=FALSE, Dbeta=TRUE, Dgamma=TRUE, Dlambda=FALSE)
                ## ev.phi = svd(-obj.cur$Hes.phi)$d
                ## ev.beta = svd(-obj.cur$Hes.beta)$d
                ## levidence = obj.cur$lpen -.5*sum(log(ev.beta[ev.beta>1-6]))
                levidence = obj.cur$lpen -.5*ldet.fun(-obj.cur$Hes.beta)
                if (!nogamma){
                    ## ev.gamma = svd(-obj.cur$Hes.gamma)$d
                    ## levidence = levidence -.5*sum(log(ev.gamma[ev.gamma>1-6]))
                    levidence = levidence -.5*ldet.fun(-obj.cur$Hes.gamma)
                }
                ## levidence = obj.cur$lpen -.5*sum(log(ev.phi[ev.phi>1-6])) -.5*sum(log(ev.beta[ev.beta>1-6])) -.5*sum(log(ev.gamma[ev.gamma>1-6]))
                ## ev.regr = svd(-obj.cur$Hes.regr)$d
                ## levidence = obj.cur$lpen -.5*sum(log(ev.phi[ev.phi>1-6])) -.5*sum(log(ev.regr[ev.regr>1-6]))
                ans = -levidence
                return(ans)
            } ## End loglambda.loss
            ##
            ## Minimize the loss function to select <log(lambda)>
            obj.ml = nlminb(start=log(c(lambda1.cur,lambda2.cur)),
                            objective=loglambda.loss,
                            lower=rep(0,J1+J2),upper=(rep(10,J1+J2)))
            lambda.cur = exp(obj.ml$par) ## Selected <lambda> --> lambda.hat
            if (J1 > 0) lambda1.cur = lambda.cur[1:J1]
            if (J2 > 0) lambda2.cur = lambda.cur[(J1+1):(J1+J2)]
        } ## Endif lambda.method == "nlminb"
        ##
        ## toc() ; cat("\n")
        ##
        ## ======================================
        ## -4- Status at the end of the iteration
        ## ======================================
        ##
        obj.cur = ff(phi.cur, beta.cur, gamma.cur,
                  tau=tau.cur, lambda1=lambda1.cur, lambda2=lambda2.cur,
                  Dphi=TRUE, D2phi=TRUE, Dbeta=TRUE, Dgamma=!nogamma)
        ##
        if (nogamma) obj.cur$grad.gamma = 0
        grad.cur = with(obj.cur, grad.regr) ##c(grad.beta,grad.gamma))
        grad.L2 = with(obj.cur, c(beta=L2norm(abs(grad.beta)),gamma=L2norm(abs(grad.gamma))))
        ##
        lpen.final = obj.cur$lpen ## Value of the penalized llik at the end of the global iteration
        ## ctrl = abs(lpen.final-lpen.start) / (abs(lpen.final) + .1)
        ## converged = (ctrl < 1e-5)
        ## converged = all(abs(grad.cur) < grad.tol*10)
        ##
        ## Effective dims & Information criteria
        ED.cur = ED.fun(obj.cur,Wood.test=Wood.test)
        ED1.tot = regr1$nfixed + ifelse(J1==0, 0, sum(ED.cur$ED1[,1]))
        ED2.tot = ifelse(nogamma, 0, regr2$nfixed) + ifelse(J2==0, 0, sum(ED.cur$ED2[,1]))
        ED.tot = ED1.tot + ED2.tot
        AIC = obj.cur$dev + 2*ED.tot
        BIC = obj.cur$dev + log(sum(n.event))*ED.tot
        ## log(evidence)
        ev.phi = svd(-obj.cur$Hes.phi)$d
        ev.beta = svd(-obj.cur$Hes.beta)$d
        if (!nogamma){
            ev.gamma = svd(-obj.cur$Hes.gamma)$d
        } else {
            ev.gamma = 0
        }
        levidence = obj.cur$lpen -.5*sum(log(ev.beta[ev.beta>1-6]))
        if (!nogamma) levidence = levidence -.5*sum(log(ev.gamma[ev.gamma>1-6]))
        ## ev.regr = svd(-obj.cur$Hes.regr)$d
        ## levidence = obj.cur$lpen -.5*sum(log(ev.phi[ev.phi>1-6])) -.5*sum(log(ev.regr[ev.regr>1-6]))
        ##
        if (iter.verbose){
            ## cat(iter,": Dev:",round(obj.cur$dev,2),
            ##     "; lpen:",obj.cur$lpen," ; AIC:",AIC," ; BIC:",BIC,
            ##     "; levidence:",levidence,"\n")
            cat(iter,
                ": levidence:",round(levidence,2),
                "; BIC:",round(BIC,2),
                "; AIC:",round(AIC,2),
                "; Dev:",round(obj.cur$dev,2),
                "; lpen:",round(obj.cur$lpen,2),
                "\n")
            ## cat("lambda1:",lambda1.cur," ; ","lambda2:",lambda2.cur,"\n")
            ## cat("ED1:",ED1.tot," ; ","ED2:",ED2.tot,"\n")
            ## cat("ED:") ; print(ED.cur)
        }
        ##
        switch(criterion, ## CONVERGENCE criterion
               "levidence" = {
                   converged = (abs(levidence-levidence.old) < criterion.tol)
                   levidence.old = levidence
               },
               "deviance" = {
                   dev.cur = obj.cur$dev
                   converged = (abs(dev.old-dev.cur) < criterion.tol)
                   dev.old = dev.cur
               },
               "lpen" = {
                   converged = (abs(lpen.old-obj.cur$lpen) < criterion.tol)
                   lpen.old = obj.cur$lpen
               },
               "AIC" = {
                   converged = (abs(AIC.old-AIC) < criterion.tol)
                   AIC.old = AIC
               },
               "BIC" = {
                   converged = (abs(BIC.old-BIC) < criterion.tol)
                   BIC.old = BIC
               },
               "gradient" = {
                   converged = all(grad.L2 < 10*grad.tol) & (iter > 5)
               }
               )
        if (iter >= iterlim) break
    } ## End While (global estimation loop)
    ##
    ## Prepare final output (after convergence)
    ## ========================================
    fit = obj.cur
    fit$criterion = criterion
    ## Report parameter estimates with their se's, z-score, Pval
    fun = function(est,se){
        mat = cbind(est=est,se=se,
                    low=est-z.alpha*se, up=est+z.alpha*se,
                    "Z"=est/se,
                    "Pval"=1-pchisq((est/se)^2,1))
        attr(mat,"ci.level") = ci.level
        return(mat)
    } ## End fun
    ##
    fit$grad.phi = obj.cur$grad.phi
    fit$grad.psi = obj.cur$grad.phi[-k.ref]
    fit$Hes.phi0 = obj.cur$Hes.phi0 ; fit$Hes.phi = obj.cur$Hes.phi
    fit$Hes.psi0 = obj.cur$Hes.phi0[-k.ref,-k.ref] ; fit$Hes.psi = obj.cur$Hes.phi[-k.ref,-k.ref]
    se.psi = sqrt(diag(solve(-fit$Hes.psi)))
    se.phi = psi2phi(se.psi)
    ##
    se.beta = sqrt(diag(solve(-fit$Hes.beta)))
    if (nogamma){
        se.gamma = 0
    } else {
        se.gamma = sqrt(diag(solve(-fit$Hes.gamma)))
    }
    ## The following line of code was there before, yielding misleading results !!
    ## se.gamma = ifelse(nogamma, 0, sqrt(diag(solve(-fit$Hes.gamma))))
    ##
    fit$gam = fit$gamma
    fit$phi = with(fit, fun(est=phi,se=se.phi))
    fit$tau = tau.cur ; fit$pen.order0 = pen.order0
    fit$beta = with(fit, fun(est=beta,se=se.beta))
    fit$gamma = with(fit, fun(est=gamma,se=se.gamma))
    if (J1 > 0){
        fit$lambda1 = lambda1.cur
        if ((lambda.method == "Schall") || (lambda.method == "LPS2")){
            fit$xi1 = log(lambda1.cur-lambda1.min)
            fit$Hes.xi1 = Hes.xi1
        }
        fit$pen.order1 = pen.order1
    }
    if (J2 > 0){
        fit$lambda2 = lambda2.cur
        if ((lambda.method == "Schall") || (lambda.method == "LPS2")){
            fit$xi2 = log(lambda2.cur-lambda2.min)
            fit$Hes.xi2 = Hes.xi2
        }
        fit$pen.order2 = pen.order2
    }
    fit$tau.method = tau.method
    fit$lambda.method = lambda.method
    ##
    ## Evaluate EDF, Chi2 and Pval for estimated additive terms
    temp = ED.fun(fit,Wood.test=Wood.test)
    fit$ED1 = temp$ED1 ; fit$ED2 = temp$ED2
    fit$ED1.Tr = temp$ED1.Tr ; fit$ED2.Tr = temp$ED2.Tr
    fit$ED1.Chi2 = temp$ED1.Chi2 ; fit$ED2.Chi2 = temp$ED2.Chi2
    ##
    ## Effective total number of parameters in regression submodels
    ED1.tot = ED2.tot = 0
    ED1.tot = regr1$nfixed + ifelse(J1==0, 0, sum(fit$ED1[,1]))
    ED2.tot = ifelse(nogamma, 0, regr2$nfixed) + ifelse(J2==0, 0, sum(fit$ED2[,1]))
    ED.tot = ED1.tot + ED2.tot
    ##
    ## AIC and BIC
    fit$nobs = nrow(regr1$Xcal)
    fit$n = n ## Numer of units (not to be confused with the number <nobs> of longitudinal binary observations)
    fit$d = sum(n.event) ## Total number of events observed on the <n> units
    ##
    fit$ED1.tot = ED1.tot ; fit$ED2.tot = ED2.tot ; fit$ED.tot = ED.tot
    fit$AIC = fit$dev + 2*ED.tot
    fit$BIC = fit$dev + log(fit$d)*ED.tot
    ## fit$BIC = fit$dev + log(fit$nobs)*ED.tot
    ##
    ## ## levidence
    ev.phi = svd(-fit$Hes.phi)$d
    ev.regr = svd(-fit$Hes.regr)$d
    fit$levidence = fit$lpen -.5*sum(log(ev.phi[ev.phi>1-6])) -.5*sum(log(ev.regr[ev.regr>1-6]))
    ##
    fit$iter = iter
    fit$elapsed.time <- (proc.time()-ptm)[1] ## Elapsed time
    ##
    ans = list(formula1=formula1, formula2=formula2, method=method,
            regr1=regr1,regr2=regr2, K0=K0,
            fit=fit,
            call=cl,
            converged=converged
            )
    ##
    ans$logLik = ans$fit$llik
    class(ans) = "tvcure"
    return(ans)
}
## End tvcure
