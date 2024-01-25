#' Predict method for tvcure model fits
#' @description
#' Predicted values based on a tvcure object.
#'
#' @usage \method{predict}{tvcure}(x, df.new, ci.level=.95, ...)
#'
#' @param x A \code{\link{tvcure.object}}.
#' @param df.new A data frame in which to look for the 'id' (distinguishing the different units), 'time' and covariate values for which 'predictions' should be made. Time values for a given 'id' should be a series of consecutive integers starting with 1. If \code{df.new$id} does not exist, then predictions are assumed to concern a single unit with consecutive time values starting with 1.
#' @param ci.level Credible level for the reported estimates. (Default: 0.95).
#' @param ... additional generic arguments.
#'
#'
#' @return A list containing, in addition to the optional \code{df.new} entries, the following elements:
#' \itemize{
#' \item{\code{Hp} : \verb{ }}{Matrix containing estimates of the cumulative population hazard \eqn{H_p(t|x_{1:t})} with its credible interval bounds at time \eqn{t} given the history of covariates.}
#' \item{\code{lHp} : \verb{ }}{Matrix containing estimates of the log cumulative population hazard \eqn{\log H_p(t|x_{1:t})} with its standard error and credible interval bounds at time \eqn{t} given the history of covariates.}
#' \item{\code{se.lHp} : \verb{ }}{Vector containing the standard errors of the estimated log cumulative population hazard at time \eqn{t} given the history of covariates.}
#' \item{\code{hp} : \verb{ }}{Matrix containing estimates of the population hazard \eqn{h_p(t|x_{1:t})} with its credible interval bounds at time \eqn{t} given the history of covariates.}
#' \item{\code{lhp} : \verb{ }}{Matrix containing estimates of the log population hazard \eqn{\log h_p(t|x_{1:t})} with its standard error and credible interval bounds at time \eqn{t} given the history of covariates.}
#' \item{\code{se.lhp} : \verb{ }}{Vector containing the standard errors of the estimated log population hazard at time \eqn{t} given the history of covariates.}
#' \item{\code{Sp} : \verb{ }}{Matrix containing estimates of the population survival fuction \eqn{S_p(t|x_{1:t})=\exp(-H_p(t|x_{1:t}))} with its credible interval bounds at time \eqn{t} given the history of covariates.}
#' \item{\code{pcure} : \verb{ }}{Matrix containing estimates of the conditional cure probability of a unit still at tisk at time \eqn{t}, \eqn{P(T=+\infty|T>t,x=x_t)}, with its credible interval bounds at time \eqn{t} if covariates remain constant from time \eqn{t}.}
#' \item{\code{llpcure} : \verb{ }}{Matrix containing estimates of the conditional log-log cure probability of a unit still at tisk at time \eqn{t}, \eqn{\log(-\log P(T=+\infty|T>t,x=x_t))}, with its standard error and credible interval bounds at time \eqn{t} if covariates remain constant from time \eqn{t}.}
#' \item{\code{se.llpcure} : \verb{ }}{Vector containing the standard errors of the estimated conditional log-log cure probability of a unit still at tisk at time \eqn{t}, \eqn{\log(-\log P(T=+\infty|T>t,x=x_t))}, if covariates remain constant from time \eqn{t}.}
#' }
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#' @references Lambert, P. and Kreyenfeld, M. (2024). Exogenous time-varying covariates in double additive cure survival model
#' with application to fertility. \emph{Journal of the Royal Statistical Society, Series A}, under review.
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
#' tau.0 = 2.7 ; lambda1.0 = c(40,15) ; lambda2.0 = c(25,70) ## Optional
#' model = tvcure(~z1+z2+s(x1)+s(x2), ~z3+z4+s(x3)+s(x4), df=df.raw,
#'                tau.0=tau.0, lambda1.0=lambda1.0, lambda2.0=lambda2.0)
#'
#' ## Covariate profiles for which 'predicted' values are requested
#' df.new = subset(df.raw, id==1 | id==4)[,-3] ## Focus on units 1 & 4
#' pred = predict(model,df.new)
#'
#' ## Visualize the estimated population survival fns for units 1 & 4
#' ## par(mfrow=c(1,2))
#' with(subset(pred,id==1), plotRegion(time,Sp,main="Id=1",
#'                               ylim=c(0,1),xlab="t",ylab="Sp(t)"))
#' with(subset(pred,id==4), plotRegion(time,Sp,main="Id=4",
#'                               ylim=c(0,1),xlab="t",ylab="Sp(t)"))
predict.tvcure <- function(x, df.new, ci.level=.95, ...){
    obj.tvcure = x
    ## Check that <id> entry in df.new. If missing, create one
    if (is.null(df.new$id)) df.new$id = rep(1,nrow(df.new))
    ##
    obj = obj.tvcure
    baseline = obj$baseline
    t.grid = obj$fit$t.grid
    ##
    f0.grid = obj$fit$f0.grid
    F0.grid = obj$fit$F0.grid
    S0.grid = obj$fit$S0.grid
    ##
    dlf0.grid = obj$fit$dlf0.grid
    dlF0.grid = obj$fit$dlF0.grid
    dlS0.grid = obj$fit$dlS0.grid
    ##
    id = df.new$id
    ids = unique(id)
    ans = df.new
    ans$Sp = ans$lhp = ans$se.lhp = ans$Hp = ans$se.lHp = numeric(nrow(df.new))
    ans$llpcure = ans$pcure = ans$se.llpcure = numeric(nrow(df.new))
    for (ii in 1:length(ids)){ ## Loop over df.new$id values
        df.sub = subset(df.new, id==ids[ii])
        idx = which(id==ids[ii])
        time = df.sub$time
        if (!all(time %in% t.grid) | min(time)!=1 | !all(diff(time) == 1)){
            cat("<df.new$time> for a given <id> should be a sequence of consecutive integers in ",min(t.grid),":",max(t.grid),"\n",sep="")
            if (min(time)!=1) cat("In addition, min(df.new$time) for a given <id> must be equal to 1\n")
            return(NULL)
        }
        ## Estimated regression and spline parameters
        phi = obj$fit$phi
        beta = obj$fit$beta
        gamma = obj$fit$gamma
        ## Design matrices for the new data
        regr1 = DesignFormula(obj$formula1, data=df.sub, K=obj$regr1$K, pen.order=obj$regr1$pen.order, knots.x=obj$regr1$knots.x)
        regr2 = DesignFormula(obj$formula2, data=df.sub, K=obj$regr2$K, pen.order=obj$regr2$pen.order, knots.x=obj$regr2$knots.x, nointercept=TRUE)
        ## Check if <gamma> is not simply constant (and equal to 1)
        q2 = ncol(regr2$Xcal) ## Total number of regression and spline parameters in short-term survival
        if (is.null(q2)){
            constant.gamma = TRUE ## Check whether covariates are absent in <formula2> --> gamma=0
        } else {
            constant.gamma = FALSE
        }
        ## Linear predictors for the new data
        eta.1 = c(regr1$Xcal %*% beta[,"est"])  ## Linear predictors in long-term survival
        if (constant.gamma){ ## Absence of covariate (and no intercept for identification reasons)
            eta.2 = 0
        } else {
            eta.2 = c(regr2$Xcal %*% gamma[,"est"]) ## Linear predictors in short-term survival
        }
        ## Fitted log(h_p(t|x)) and log(H_p(t|x))
        ## --------------------------------------
        Delta = 1
        switch(baseline,
               "F0" = {
                   lhp = eta.1 + eta.2 + (exp(eta.2)-1) * log(F0.grid[time]) + log(f0.grid[time])
                   Hp = cumsum(exp(lHp)*Delta)
                   lHp = log(Hp)
                   ## lHp = eta.1 + exp(eta.2)*log(F0.grid[time])
               },
               "S0" = {
                   lhp = eta.1 + eta.2 + (exp(eta.2)-1) * log(S0.grid[time]) + log(f0.grid[time]) ## <-- Pack_v2
                   Hp = cumsum(exp(lhp)*Delta)
                   lHp = log(Hp)
                   ## lHp = eta.1 + log(1-(S0.grid[time]^exp(eta.2))) ## <-- Pack_v2
               }
               )
        ## Computation of se(log(Hp)) & se(log(hp))
        ## ----------------------------------------
        q2 = ncol(regr2$Xcal) ## Total number of regression and spline parameters in short-term survival
        if (is.null(q2)){
            constant.gamma = TRUE ## Check whether covariates are absent in <formula2> --> gamma=0
        } else {
            constant.gamma = FALSE
        }
        ## Derivatives of log(Hp) & log(hp) wrt <beta> at the grid times
        Dlhp.beta = regr1$Xcal ## Tgrid x nbeta
        DlHp.beta = (1/Hp) * apply(exp(lhp) * regr1$Xcal * Delta, 2, cumsum) ## Tgrid x nbeta
        dim(DlHp.beta) = dim(regr1$Xcal) ## ... to force matrix format if single row matrix
        ##
        ## Derivatives of log(Hp) & log(hp) wrt <gamma> at the grid times
        if (!constant.gamma){
            switch(baseline,
                   "F0" = {
                       tempo = exp(eta.2)*log(F0.grid[time])
                   },
                   "S0" = {
                       tempo = exp(eta.2)*log(S0.grid[time]) ## <-- Pack_v2
                   }
                   )
            Dlhp.gamma = (1+tempo) * regr2$Xcal ## Tgrid x ngamma
            DlHp.gamma = (1/Hp) * apply(exp(lhp) * (1+tempo) * regr2$Xcal * Delta, 2, cumsum) ## Tgrid x ngamma
            dim(DlHp.gamma) = dim(regr2$Xcal) ## ... to force matrix format if single row matrix

        }
        ## Derivatives of log(Hp) & log(hp) wrt <beta, gamma> at the grid times
        if (!constant.gamma){
            Dlhp.regr = cbind(Dlhp.beta, Dlhp.gamma) ## Tgrid x (nbeta + ngamma)
            DlHp.regr = cbind(DlHp.beta, DlHp.gamma) ## Tgrid x (nbeta + ngamma)
        } else {
            Dlhp.regr = Dlhp.beta
            DlHp.regr = DlHp.beta
        }
        ##
        ## Derivatives of log(Hp) & log(hp) wrt <psi> at the grid times with <psi> = <phi[-k.ref]>
        k.ref = obj$fit$k.ref
        npsi = nrow(obj$fit$phi) - 1
        if (!constant.gamma){
            switch(baseline,
                   "F0" = {
                       part1 = (exp(eta.2)-1) * dlF0.grid[time,-k.ref,drop=FALSE]
                   },
                   "S0" = {
                       part1 = (exp(eta.2)-1) * dlS0.grid[time,-k.ref,drop=FALSE]
                   }
           )
        } else {
            part1 = matrix(0, nrow=length(time), ncol=npsi)
        }
        part2 = dlf0.grid[time,-k.ref,drop=FALSE]
        ##
        Dlhp.psi = part1 + part2  ## Tgrid x (K0-1)
        DlHp.psi = (1/Hp) * apply(exp(lhp) * Dlhp.psi * Delta, 2, cumsum) ## Tgrid x (K0-1)
        dim(DlHp.psi) = dim(Dlhp.psi) ## To force a matrix format after "apply" (even with single row)
        ##
        ## se(log(hp)) wrt <beta, gamma, psi>
        Diag.1 = diag(Dlhp.regr %*% solve(-obj$fit$Hes.regr) %*% t(Dlhp.regr))
        Diag.1 = Diag.1 + diag(Dlhp.psi %*% solve(-obj$fit$Hes.psi) %*% t(Dlhp.psi))
        se.lhp = sqrt(Diag.1)
        ##
        ## se(log(Hp)) wrt <beta, gamma, psi>
        Diag.2 = diag(DlHp.regr %*% solve(-obj$fit$Hes.regr) %*% t(DlHp.regr))
        Diag.2 = Diag.2 + diag(DlHp.psi %*% solve(-obj$fit$Hes.psi) %*% t(DlHp.psi))
        se.lHp = sqrt(Diag.2)
        ##
        ## Cure probability pcure(t|x) = P(T=inf | T>t & X=x after t)
        ##    log(-log pcure(t|x)): llpcure
        ##    pcure(t|x) = exp(-exp(llpcure))
        ## ----------------------------------------------------------
        switch(baseline,
               "F0" = {
                   llpcure = eta.1 + log(1-F0.grid[time]^(exp(eta.2)))
               },
               "S0" = {
                   llpcure = eta.1 + exp(eta.2) * log(S0.grid[time])
               }
        )
        ##
        ## se(llpcure)
        if (baseline == "S0"){
            Dllpcure.beta = regr1$Xcal ## Tgrid x nbeta
            if (!constant.gamma){
                Dllpcure.gamma = tempo * regr2$Xcal ## Tgrid x ngamma
                Dllpcure.regr = cbind(Dllpcure.beta, Dllpcure.gamma)
            } else {
                Dllpcure.regr = Dllpcure.beta
            }
            Dllpcure.psi = exp(eta.2) * dlS0.grid[time,-k.ref,drop=FALSE]
        }
        Diag.3 = diag(Dllpcure.regr %*% solve(-obj$fit$Hes.regr) %*% t(Dllpcure.regr))
        Diag.3 = Diag.3 + diag(Dllpcure.psi %*% solve(-obj$fit$Hes.psi) %*% t(Dllpcure.psi))
        se.llpcure = sqrt(Diag.3)
        ##
        ## Returned quantities per unit
        ## ----------------------------
        ans$llpcure[idx] = llpcure ## log(-log P(cure))
        ans$se.llpcure[idx] = se.llpcure
        ## ans$pcure[idx] = exp(-exp(llpcure)) ## P(cure)
        ans$lhp[idx] = lhp ## log hp(t)
        ans$se.lhp[idx] = se.lhp
        ans$Hp[idx] = Hp ## Hp(t)
        ans$se.lHp[idx] = se.lHp
        ans$Sp[idx] = exp(-Hp) ## Sp(t)
        ## ans$Su[idx] = (exp(-Hp) - exp(-exp(eta.1))) / (1-exp(-exp(eta.1)))
    }
    ## Final returned objects
    ## ----------------------
    z = qnorm(1-.5*(1-ci.level))
    ans$llpcure = with(ans, cbind(est=llpcure, se=se.llpcure, low=llpcure-z*se.llpcure, up=llpcure+z*se.llpcure))
    ans$pcure = exp(-exp(ans$llpcure[,c(1,4,3),drop=FALSE])) ; colnames(ans$pcure) = c("est","low","up")
    ans$lhp = with(ans, cbind(est=lhp, se=se.lhp, low=lhp-z*se.lhp, up=lhp+z*se.lhp))
    ans$hp = exp(ans$lhp[,c(1,3,4),drop=FALSE])
    ans$lHp = with(ans, cbind(est=log(Hp), se=se.lHp, low=log(Hp)-z*se.lHp, up=log(Hp)+z*se.lHp))
    ans$Hp = exp(ans$lHp[,c(1,3,4),drop=FALSE])
    ans$Sp = exp(-ans$Hp[,c(1,3,2)]) ; colnames(ans$Sp) = colnames(ans$Hp)[c(1,3,2)]
    ##
    return(ans)
}



## predictOld.tvcure <- Function(obj.tvcure, df.new){
##     ## Check that <id> entry in df.new. If missing, create one
##     if (is.null(df.new$id)) df.new$id = rep(1,nrow(df.new))
##     ##
##     obj = obj.tvcure
##     baseline = obj$baseline
##     t.grid = obj$fit$t.grid
##     f0.grid = obj$fit$f0.grid
##     F0.grid = obj$fit$F0.grid
##     S0.grid = obj$fit$S0.grid
##     ##
##     id = df.new$id
##     ids = unique(id)
##     ans = df.new
##     ans$Sp = ans$Su = ans$Hp = ans$lhp = rep(NA,nrow(df.new))
##     for (ii in 1:length(ids)){ ## Loop over df.new$id values
##         df.sub = subset(df.new, id==ids[ii])
##         idx = which(id==ids[ii])
##         time = df.sub$time
##         if (!all(time %in% t.grid) | min(time)!=1 | !all(diff(time) == 1)){
##             cat("<df.new$time> for a given <id> should be a sequence of consecutive integers in ",min(t.grid),":",max(t.grid),"\n",sep="")
##             if (min(time)!=1) cat("In addition, min(df.new$time) for a given <id> must be equal to 1\n")
##             return(NULL)
##         }
##         ## Estimated regression and spline parameters
##         phi = obj$fit$phi
##         beta = obj$fit$beta
##         gamma = obj$fit$gamma
##         ## Design matrices for the new data
##         regr1 = DesignFormula(obj$formula1, data=df.sub, K=obj$regr1$K, pen.order=obj$regr1$pen.order, knots.x=obj$regr1$knots.x)
##         regr2 = DesignFormula(obj$formula2, data=df.sub, K=obj$regr2$K, pen.order=obj$regr2$pen.order, knots.x=obj$regr2$knots.x)
##         ## Check if <gamma> is not simply constant (and equal to 1)
##         q2 = ncol(regr2$Xcal) ## Total number of regression and spline parameters in short-term survival
##         if (is.null(q2)){
##             constant.gamma = TRUE ## Check whether covariates are absent in <formula2> --> gamma=0
##         } else {
##             constant.gamma = FALSE
##         }
##         ## Linear predictors for the new data
##         eta.1 = c(regr1$Xcal %*% beta[,"est"])  ## Linear predictors in long-term survival
##         if (constant.gamma){ ## Absence of covariate (and no intercept for identification reasons)
##             eta.2 = 0
##         } else {
##             eta.2 = c(regr2$Xcal %*% gamma[,"est"]) ## Linear predictors in short-term survival
##         }
##         ## Fitted log(h_p(t|x)) and log(H_p(t|x))
##         Delta = 1
##         switch(baseline,
##                "F0" = {
##                    lhp = eta.1 + eta.2 + (exp(eta.2)-1) * log(F0.grid[time]) + log(f0.grid[time])
##                    Hp = cumsum(exp(lHp)*Delta)
##                    ## lHp = eta.1 + exp(eta.2)*log(F0.grid[time])
##                },
##                "S0" = {
##                    lhp = eta.1 + eta.2 + (exp(eta.2)-1) * log(S0.grid[time]) + log(f0.grid[time]) ## <-- Pack_v2
##                    Hp = cumsum(exp(lhp)*Delta)
##                    ## lHp = eta.1 + log(1-(S0.grid[time]^exp(eta.2))) ## <-- Pack_v2
##                }
##                )
##         ##
##         ans$Su[idx] = (exp(-Hp) - exp(-exp(eta.1))) / (1-exp(-exp(eta.1)))
##         ##
##         ans$lhp[idx] = lhp
##         ans$Hp[idx] = Hp
##         ans$Sp[idx] = exp(-Hp)
##     }
##     ##
##     return(ans)
## }
