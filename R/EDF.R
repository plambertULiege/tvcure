#' Compute the effective degrees freedom in a tvcure model
#'
#' @param model A tvcure object
#' @param Wood.test Logical indicating if P-values based on Wood's test (Biometrika 2013) of the significance of additive terms should be preferred over basic Chi-square tests. (Default: FALSE).
#' @param joint.computation Logical indicating if variance-covariance matrices for the regression and spline parameters in the long- and short-term survival submodels should be computed jointly (TRUE) or separately (FALSE). (Default: TRUE).
#'
#' @return A list containing the effective degrees of freedom for the additive terms in the long-term (quantum) and short-term (timing) survival submodels, with the selected statistical test for significance and its P-value.
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#' @references Lambert, P. and Kreyenfeld, M. (2025).
#' Time-varying exogenous covariates with frequently changing values in double additive cure survival model: an application to fertility.
#' \emph{Journal of the Royal Statistical Society, Series A}. <doi:10.1093/jrsssa/qnaf035>
#'
#' @examples
#' \donttest{
#' require(tvcure)
#' ## Simulated data generation
#' beta = c(beta0=.4, beta1=-.2, beta2=.15) ; gam = c(gam1=.2, gam2=.2)
#' data = simulateTVcureData(n=500, seed=123, beta=beta, gam=gam,
#'                           RC.dist="exponential",mu.cens=550)$rawdata
#' ## TVcure model fitting
#' tau.0 = 2.7 ; lambda1.0 = c(40,15) ; lambda2.0 = c(25,70) ## Optional
#' model = tvcure(~z1+z2+s(x1)+s(x2), ~z3+z4+s(x3)+s(x4), data=data,
#'                tau.0=tau.0, lambda1.0=lambda1.0, lambda2.0=lambda2.0)
#' EDF(model)
#' }
#' @export
EDF = function(model, Wood.test=FALSE, joint.computation=TRUE){
    fit = model$fit
    ## marginalized = fit$marginalised ; if (is.null(marginalized)) marginalized = FALSE
    marginalized = ifelse(exists("marginalized",fit), fit$marginalized, FALSE)
    Tr.test = function(p,Xt,V,edf){
        ans <- tvcure::testStat(p, Xt, V, min(ncol(Xt), edf), type = 0, res.df = -1)
        return(ans)
    }
    ED1 = ED2 = ED1.Tr = ED2.Tr = ED1.Chi2 = ED2.Chi2 = ED1.linear = ED2.linear = NULL
    ##
    nbeta = fit$nbeta   ## Nbr of regression & spline parameters in long-term survival
    ngamma = fit$ngamma ## Nbr of regression & spline parameters in short-term survival
    nogamma = fit$nogamma ## Indicates when covariates are absent in <formula2>
    J1 = fit$J1 ; J2 = fit$J2 ; K1 = fit$K1 ; K2 = fit$K2
    nfixed1 = fit$nfixed1 ; nfixed2 = fit$nfixed2
    ##
    Sigma.regr = solve(-fit$Hes.regr+1e-6*diag(ncol(fit$Hes.regr)))
    Sigma.regr0 = solve(-fit$Hes.regr0+1e-6*diag(ncol(fit$Hes.regr)))
    if (joint.computation){
        ED.full = rowSums(t(Sigma.regr) * (-fit$Hes.regr0))
        idx1 = 1:nbeta
        Sigma.beta = Sigma.regr[idx1,idx1]
        ED.beta = ED.full[idx1]
    } else {
        Sigma.beta = with(fit, solve(-Hes.beta+1e-6*diag(ncol(Hes.beta))))
        ED.beta = rowSums(t(Sigma.beta) * (-fit$Hes.beta0))
    }
    ##
    if (!nogamma){
        if (joint.computation){
            Sigma.gamma = Sigma.regr[-idx1,-idx1]
            ED.gamma = ED.full[-idx1]
        } else {
            Sigma.gamma = with(fit, solve(-Hes.gamma+1e-6*diag(ncol(Hes.gamma))))
            ED.gamma = rowSums(t(Sigma.gamma) * (-fit$Hes.gamma0))
        }
    } else {
        Sigma.gamma = NULL
        ED.gamma = 0
    }
    ##
    ED.full = c(ED.beta,ED.gamma) ## Added on 2023.10.11
    ##
    ## Sigma.regr = with(fit, solve(-Hes.regr)) ## Removed on 2023.10.11
    ## ED.full = rowSums(t(Sigma.regr) * (-fit$Hes.regr0)) ## Removed on 2023.10.11
    ##
    ED1.CI = ED2.CI = NULL
    if (J1 > 0){
        ## Sigma.beta = with(fit, solve(-Hes.beta))
        ## ED1.full = with(fit, rowSums(t(Sigma.beta) * (-Hes.beta0)))
        ## ED1.full = with(fit, diag(Sigma.beta %*% (-Hes.beta0)))
        ED1 = Chi2 = Tr = Pval.Tr = Pval.Chi2 = Chi2.lin = Pval.lin = numeric(J1)
        ED1.CI = matrix(nrow=J1,ncol=4) ; rownames(ED1.CI) = paste("f1(",fit$additive.lab1,")",sep="")
        ## colnames(ED1.CI) = c("edf","low","up","Plin")
        colnames(ED1.CI) = c("edf","low","up","P(<1.5)")
        for (j in 1:J1){
            idx = nfixed1 + (j-1)*K1 + (1:K1)
            beta.j = fit$beta[idx]
            if (joint.computation){
                ED1[j] = sum(ED.full[idx])
            } else {
                Sigma.betaj = with(fit, solve(-Hes.beta[idx,idx]+1e-6*diag(length(idx))))
                bele = t(Sigma.betaj) * (-fit$Hes.beta0[idx,idx])
                ED1[j] = sum(bele)
            }
            if (marginalized) ED1[j] = model$fit$margins1[[j]]$mu.edf ## Posterior mean
            ## Wood's significance test
            ## ------------------------
            ngrid = 200
            knots.x = fit$knots1.x[[j]] ## knots.x = regr1$knots.x[[j]]
            xL = min(knots.x) ; xU = max(knots.x)
            pen.order = fit$pen.order1 ## regr1$pen.order
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
            ## Chi2 significance test
            ## ----------------------
            Chi2[j] = sum(beta.j * c(solve(Sigma.beta[idx,idx], beta.j))) ## Added on 2023.10.11
            ## Chi2[j] = sum(beta.j * c(solve(Sigma.regr[idx,idx]) %*% beta.j)) ## Removed on 2023.10.11
            Pval.Chi2[j] = 1 - pchisq(Chi2[j],ED1[j])
            ## cat("Naive:",Chi2[j],Pval.Chi2[j],ED1[j],"\n")
            ## End Chi2
            ##
            ## Test for linearity
            ## ------------------
            Dd = fit$Dd1.x ; theta.j = beta.j ; Sig = Sigma.beta[idx,idx] ; ED.j = ED1[j]
            ##
            Dtheta.j = c(Dd %*% theta.j) ; DSigDt = Dd %*% (Sig %*% t(Dd))
            Chi2.lin[j] = sum(Dtheta.j * solve(DSigDt, Dtheta.j))
            Pval.lin[j] = 1 - pchisq(Chi2.lin[j], ED.j-1)
            ## End test for Linearity
            ##
            ## Credible interval for ED1
            ## -------------------------
            if (marginalized){
                alpha = 1-model$fit$ci.level
                CI = model$fit$margins1[[j]]$qedf(c(.5*alpha,1-.5*alpha))
                Plin = model$fit$margins1[[j]]$pedf(1.5)
                ED1.CI[j,] = c(ED1[j],CI,Plin)
            }
        }
        ED1.Tr = cbind(ED1,Tr,Pval.Tr)
        ED1.Chi2 = cbind(ED1,Chi2,Pval.Chi2)
        ED1.linear = cbind(ED1,Chi2.lin,Pval.lin)
        rownames(ED1.Tr) = rownames(ED1.Chi2) = paste("f1(",fit$additive.lab1,")",sep="")
        colnames(ED1.Tr) = c("edf","Tr","Pval")
        colnames(ED1.Chi2) = c("edf","Chi2","Pval")
        colnames(ED1.linear) = c("edf","Chi2.lin","Pval")
    }
    ##
    if (J2 > 0){
        ## Sigma.gamma = with(fit, solve(-Hes.gamma))
        ## ED2.full = with(fit, rowSums(t(Sigma.gamma) * (-Hes.gamma0)))
        ## ED2.full = with(fit, diag(Sigma.gamma %*% (-Hes.gamma0)))
        ED2 = Chi2 = Tr = Pval.Tr = Pval.Chi2 = Chi2.lin = Pval.lin = numeric(J2)
        ED2.CI = matrix(nrow=J2,ncol=4) ; rownames(ED2.CI) = paste("f2(",fit$additive.lab2,")",sep="")
        ## colnames(ED2.CI) = c("edf","low","up","Plin")
        colnames(ED2.CI) = c("edf","low","up","P(<1.5)")
        for (j in 1:J2){
            idx = nfixed2 + (j-1)*K2 + (1:K2)
            gamma.j = fit$gamma[idx]
            if (joint.computation){
                ED2[j] = sum(ED.full[nbeta + idx])
            } else {
                Sigma.gammaj = with(fit, solve(-Hes.gamma[idx,idx]+1e-6*diag(length(idx))))
                bele = t(Sigma.gammaj) * (-fit$Hes.gamma0[idx,idx])
                ED2[j] = sum(bele)
            }
            if (marginalized) ED2[j] = model$fit$margins2[[j]]$mu.edf ## Posterior mean
            ## Wood's significance test
            ## ------------------------
            ngrid = 200
            knots.x = fit$knots2.x[[j]] ## regr2$knots.x[[j]]
            xL = min(knots.x) ; xU = max(knots.x)
            pen.order = fit$pen.order2 ## regr2$pen.order
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
            ## Chi2 significance test
            ## ----------------------
            Chi2[j] = sum(gamma.j * c(solve(Sigma.gamma[idx, idx], gamma.j))) ## Added on 2023.10.11
            ## Chi2[j] = sum(gamma.j * c(solve(Sigma.regr[nbeta+idx, nbeta+idx]) %*% gamma.j)) ## Removed on 2023.10.11
            Pval.Chi2[j] = 1 - pchisq(Chi2[j],ED2[j])
            ## cat("Naive:",Chi2[j],Pval[j],ED2[j],"\n")
            ## End Chi2
            ##
            ## Test for linearity
            ## ------------------
            Dd = fit$Dd2.x ; theta.j = gamma.j ; Sig = Sigma.gamma[idx,idx] ; ED.j = ED2[j]
            ##
            Dtheta.j = c(Dd %*% theta.j) ; DSigDt = Dd %*% (Sig %*% t(Dd))
            Chi2.lin[j] = sum(Dtheta.j * solve(DSigDt, Dtheta.j))
            Pval.lin[j] = 1 - pchisq(Chi2.lin[j], ED.j-1)
            ##
            ## Credible interval for ED2
            if (marginalized){
                alpha = 1-model$fit$ci.level
                CI = model$fit$margins2[[j]]$qedf(c(.5*alpha,1-.5*alpha))
                Plin = model$fit$margins2[[j]]$pedf(1.5)
                ED2.CI[j,] = c(ED2[j],CI,Plin)
            }
        }
        ED2.Tr = cbind(ED2,Tr,Pval.Tr)
        ED2.Chi2 = cbind(ED2,Chi2,Pval.Chi2)
        ED2.linear = cbind(ED2,Chi2.lin,Pval.lin)
        rownames(ED2.Tr) = rownames(ED2.Chi2) = paste("f2(",fit$additive.lab2,")",sep="")
        ## rownames(ED2.Tr) = rownames(ED2.Chi2) = paste("f2(",regr2$additive.lab,")",sep="")
        colnames(ED2.Tr) = c("edf","Tr","Pval")
        colnames(ED2.Chi2) = c("edf","Chi2","Pval")
        colnames(ED2.linear) = c("edf","Chi2.lin","Pval")
    }
    ##
    if ((marginalized) & (J1 > 0)){
        ED1 = cbind(ED1.CI,ED1.Tr[,-1,drop=FALSE],fit$ED1.Chi2[,-1,drop=FALSE])
    } else {
        ED1 = cbind(ED1.Tr,ED1.Chi2[,-1,drop=FALSE])
        ## ED1 = cbind(ED1.linear,ED1.Tr[,-1,drop=FALSE],ED1.Chi2[,-1,drop=FALSE])
    }
    ##
    if ((marginalized) & (J2 > 0)){
        ED2 = cbind(ED2.CI,ED2.Tr[,-1,drop=FALSE],ED2.Chi2[,-1,drop=FALSE])
    } else {
        ED2 = cbind(ED2.Tr,ED2.Chi2[,-1,drop=FALSE])
        ## ED2 = cbind(ED2.linear,ED2.Tr[,-1,drop=FALSE],ED2.Chi2[,-1,drop=FALSE])
    }
    ## if (Wood.test){
    ##     ED1=ED1.Tr ; ED2=ED2.Tr
    ## } else {
    ##     ED1=ED1.Chi2 ; ED2=ED2.Chi2
    ## }
    ##
    ED1.tot = ED2.tot = 0
    ED1.tot = nfixed1 + ifelse(J1==0, 0, sum(ED1[,1]))
    ED2.tot = ifelse(nogamma, 0, nfixed2) + ifelse(J2==0, 0, sum(ED2[,1]))
    ED.tot = ED1.tot + ED2.tot
    ##
    return(list(ED1=ED1, ED2=ED2,
                ED1.CI=ED1.CI, ED2.CI=ED2.CI,
                ED1.tot=ED1.tot, ED2.tot=ED2.tot, ED.tot=ED.tot,
                ED1.Tr=ED1.Tr, ED2.Tr=ED2.Tr,
                ED1.Chi2=ED1.Chi2, ED2.Chi2=ED2.Chi2,
                ED1.linear=ED1.linear, ED2.linear=ED2.linear
                ))
} ## End EDF
