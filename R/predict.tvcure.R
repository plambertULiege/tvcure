predict.tvcure <- function(obj.tvcure, df.new, ci.level=.95){
    ## Check that <id> entry in df.new. If missing, create one
    if (is.null(df.new$id)) df.new$id = rep(1,nrow(df.new))
    ##
    obj = obj.tvcure
    method = obj$method
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
    ans$Sp = ans$Su = ans$lhp = ans$se.lhp = ans$Hp = ans$se.lHp = numeric(nrow(df.new))
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
        switch(method,
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
            switch(method,
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
            switch(method,
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
        switch(method,
               "F0" = {
                   llpcure = eta.1 + log(1-F0.grid[time]^(exp(eta.2)))
               },
               "S0" = {
                   llpcure = eta.1 + exp(eta.2) * log(S0.grid[time])
               }
        )
        ##
        ## se(llpcure)
        if (method == "S0"){
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
    ans$Sp = exp(-ans$Hp)
    ##
    return(ans)
}



## predictOld.tvcure <- Function(obj.tvcure, df.new){
##     ## Check that <id> entry in df.new. If missing, create one
##     if (is.null(df.new$id)) df.new$id = rep(1,nrow(df.new))
##     ##
##     obj = obj.tvcure
##     method = obj$method
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
##         switch(method,
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
