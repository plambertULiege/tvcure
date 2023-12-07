#' Extract additive term estimates from a tvcure object.
#'
#' @param obj.tvcure a \code{\link{tvcure.object}}.
#' @param ngrid number of gridpoints where the fitted additive terms are evaluated.
#' @param ci.level confidence level for the pointwise credible intervals of the additive terms.
#'
#' @return A list with following elements:
#' \itemize{
#' \item{\code{f0} : \verb{ }}{a function estimate of \eqn{f_0}.}
#' \item{\code{F0} : \verb{ }}{a function estimate of \eqn{F_0}.}
#' \item{\code{T} : \verb{ }}{the follow-up time after which a unit is considered `cured'.}
#' \item{\code{nfixed1} : \verb{ }}{the number of non-penalized regression parameter in the long-term term (or quantum) submodel.}
#' \item{\code{J1} : \verb{ }}{number of additive terms in the long-term term (or quantum) submodel.}
#' \item{\code{additive.lab1} : \verb{ }}{labels of the additive terms in the long-term term (or quantum) submodel.}
#' \item{\code{K1} : \verb{ }}{number of P-spline parameters per additive term in the long-term term (or quantum) submodel.}
#' \item{\code{knots1} : \verb{ }}{list of length J1 containing the knots of the additive term in the long-term term (or quantum) submodel.}
#' \item{\code{f1.grid} : \verb{ }}{list of length J1 containing for each additive term in the long-term term (or quantum) submodel, a list of length 2 with elements <x> and <y.mat>. 
#' Element <x> is a vector of \code{ngrid} equidistant values covering the range of values for the covariate ; 
#' <y.mat> is (ngrid x 3) matrix containing in column 1 the estimated values of the additive term at <x> and the bounds of the credible interval for it in the other 2 columns.}
#' \item{\code{f1} : \verb{ }}{list of length J1 containing the estimated function of the corresponding additive term in the long-term term (or quantum) submodel.}
#' \item{\code{f1.se} : \verb{ }}{list of length J1 containing the estimated standard error function of the corresponding additive term in the long-term term (or quantum) submodel.}
#' }
#' The same definitions applies for \code{nfixed2}, \code{J2}, \code{additive.lab2}, \code{K2}, \code{knots2}, 
#' \code{f2.grid}, \code{f2}, \code{f2.se} with the additive terms in the short-term (or timing) submodel.
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
#' tau.0 = 2.5 ; lambda1.0 = c(285,15) ; lambda2.0 = c(25,1325) ## Optional
#' model = tvcure(~z1+z2+s(x1)+s(x2), ~z3+z4+s(x3)+s(x4), df=df.raw,
#'                tau.0=tau.0, lambda1.0=lambda1.0, lambda2.0=lambda2.0)
#' 
#' ## Extract additive term estimates from tvcure object
#' obj = additive.tvcure(model)
#' names(obj)
#' 
#' @export
#' 
additive.tvcure <- function(obj.tvcure,ngrid=300, ci.level=.95){
    obj = obj.tvcure
    nfixed1 = obj$regr1$nfixed ; nfixed2 = obj$regr2$nfixed ## Number of fixed parms in the regression submodels
    J1 = obj$regr1$J ; J2 = obj$regr2$J ## Number of additive terms
    K1 = obj$regr1$K ; K2 = obj$regr2$K ## Number of centered B-splines per additive term
    ## f1 = f2 = f1.se = f2.se = list() ## Fitted additive terms in long- and short-term submodels
    alpha = 1 - ci.level ; z.alpha = qnorm(1-.5*alpha)
    ans = NULL
    ## f0 & F0: dynamic for the non-cured subjects
    T = max(obj$fit$t.grid)
    f0 = function(x){
        f0 = 0*x
        idx = which((x>=0)&(x<=T))
        if (length(idx)>0) f0[idx] = splinefun(obj$fit$t.grid,obj$fit$f0.grid)(x[idx])
        return(f0)
    }
    attr(f0,"support") = c(0,T)
    ##
    F0 = function(x){
        ## ifelse(x<=0,0, ifelse(x>T,1, splinefun(obj$fit$t.grid,obj$fit$F0.grid)(x)))
        F0 = 0*x
        F0[x >= T] = 1
        idx = which((x>0)&(x<T))
        if (length(idx)>0) F0[idx] = splinefun(obj$fit$t.grid,obj$fit$F0.grid)(x[idx])
        return(F0)
    }
    attr(F0,"support") = c(0,T)
    ans$f0 = f0 ; ans$F0=F0 ; ans$T = T
    ##
    rangetol = function(x,tol=.01){
        temp = range(x) ; df = diff(temp)
        ans = c(temp[1]-.5*tol*df,temp[2]+.5*tol*df)
        return(ans)
    }
    ##
    ## Sigma.regr = with(obj$fit, solve(-Hes.regr))
    ##
    ## Additive part in long-term survival
    nbeta = obj$fit$nbeta
    ans$nfixed1=nfixed1 ; ans$J1 = J1
    if (J1 > 0){
        add.lab = obj$regr1$additive.lab
        ans$additive.lab1 = add.lab
        ans$K1=K1 ; ans$knots1 = obj$regr1$knots.x
        Sigma = with(obj$fit, solve(-Hes.beta)) ## Added in 2023.10.11
        ## Sigma = Sigma.regr[1:nbeta,1:nbeta] ## Removed in 2023.10.11
        f.grid = f = f.se = list()
        ##
        for (j in 1:J1){
            idx = nfixed1 + (j-1)*K1 + (1:K1)
            beta.j = obj$fit$beta[idx,"est"] ## Centered B-splines coefs for jth additive term
            knots.x = obj$regr1$knots.x[[j]] ; xL = min(knots.x) ; xU = max(knots.x)
            pen.order = obj$regr1$pen.order
            x.grid = seq(min(knots.x),max(knots.x),length=ngrid) ## Grid of values
            cm = obj$regr1$cm.values[[j]]
            cB = centeredBasis.gen(x.grid,knots=knots.x,cm=cm,pen.order)$B ## Centered B-spline basis
            y.grid = c(cB %*% beta.j)
            y.grid.se = sqrt(diag(cB %*% (Sigma[idx,idx]%*%t(cB))))
            ylow = y.grid - z.alpha*y.grid.se
            yup  = y.grid + z.alpha*y.grid.se
            ##
            f.grid[[add.lab[j]]]$x = x.grid
            f.grid[[add.lab[j]]]$y.mat = cbind(est=y.grid,low=ylow,up=yup)
            ##
            f[[add.lab[j]]]    = splinefun(x.grid, y.grid)
            f.se[[add.lab[j]]] = splinefun(x.grid, y.grid.se)
            attr(f[[add.lab[j]]],"support") = attr(f.se[[add.lab[j]]],"support") = c(xL,xU)
            attr(f[[add.lab[j]]],"label") = attr(f.se[[add.lab[j]]],"label") = add.lab[j]
            attr(f[[add.lab[j]]],"range") = rangetol(y.grid)
        }
        ans$f1.grid = f.grid
        ans$f1 = f ; ans$f1.se = f.se
    }
    ## Additive part in short-term survival
    ngamma = obj$fit$ngamma
    ans$nfixed2=nfixed2 ; ans$J2 = J2
    if (J2 > 0){
        add.lab = obj$regr2$additive.lab
        ans$additive.lab2 = add.lab
        ans$K2=K2 ; ans$knots2 = obj$regr2$knots.x
        Sigma = with(obj$fit, solve(-Hes.gamma)) ## Added in 2023.10.11
        ## Sigma = Sigma.regr[nbeta + (1:ngamma),nbeta + (1:ngamma)] ## Removed in 2023.10.11
        f.grid = f = f.se = list()
        for (j in 1:J2){
            idx = nfixed2 + (j-1)*K2 + (1:K2)
            gamma.j = obj$fit$gamma[idx,"est"] ## Centered B-splines coefs for jth additive term
            knots.x = obj$regr2$knots.x[[j]] ; xL = min(knots.x) ; xU = max(knots.x)
            pen.order = obj$regr2$pen.order
            x.grid = seq(min(knots.x),max(knots.x),length=ngrid) ## Grid of values
            cm = obj$regr2$cm.values[[j]]
            cB = centeredBasis.gen(x.grid,knots=knots.x,cm=cm,pen.order)$B ## Centered B-spline basis
            y.grid = c(cB %*% gamma.j)
            y.grid.se = sqrt(diag(cB %*% (Sigma[idx,idx]%*%t(cB))))
            ylow = y.grid - z.alpha*y.grid.se
            yup  = y.grid + z.alpha*y.grid.se
            ##
            f.grid[[add.lab[j]]]$x = x.grid
            f.grid[[add.lab[j]]]$y.mat = cbind(est=y.grid,low=ylow,up=yup)
            ##
            f[[add.lab[j]]]    = splinefun(x.grid, y.grid)
            f.se[[add.lab[j]]] = splinefun(x.grid, y.grid.se)
            attr(f[[add.lab[j]]],"support") = attr(f.se[[add.lab[j]]],"support") = c(xL,xU)
            attr(f[[add.lab[j]]],"label") = attr(f.se[[add.lab[j]]],"label") = add.lab[j]
            attr(f[[add.lab[j]]],"range") = rangetol(y.grid)
        }
        ans$f2.grid = f.grid
        ans$f2 = f ; ans$f2.se = f.se
    }
    ##
    return(ans)
}

# ## Generic function to plot credible region
# plotRegion = function(x,mat,add=FALSE,xlim=range(x),ylim=range(mat),
#                       lwd=2,xlab="",ylab="",main="",...){
#   f = mat[,1] ; f.low = mat[,2] ; f.up = mat[,3]
#   if (add==FALSE) plot(x,f,type="n",ylim=ylim,xlim=xlim,
#                        lwd=lwd,xlab=xlab,ylab=ylab,main=main,...)
#   polygon(c(x,rev(x)),c(f.low,rev(f.up)),col="grey",border=F)
#   lines(x,f,lwd=lwd)
# }



# ##
# levidence <- function(object, ...) UseMethod("levidence")
# levidence.tvcure <- function(obj.tvcure, ..., k=2){
#     obj = obj.tvcure
#     lls = function(obj) return(ans = c(levidence=obj$fit$levidence, edf=obj$fit$ED.tot, nobs=obj$fit$nobs))
#     if (!missing(...)) {
#         vals = sapply(list(obj,...), lls)
#         val <- data.frame(edf = round(vals[2L, ],2), levidence = vals[1L, ])
#         nos <- na.omit(vals[3L, ])
#         if (length(nos) && any(nos != nos[1L])) warning("models are not all fitted to the same number of observations")
#         Call <- match.call()
#         Call$k <- NULL
#         row.names(val) <- as.character(Call[-1L])
#         val
#     } else {
#         vals = unname(lls(obj))
#         vals[1L] + k * vals[2L]
#     }
# }

## --------------------------------
## Previous GGPLOT2 implementation
## --------------------------------
## ##
## tvcure_theme <- function(h_just = 0.5){
##     ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = h_just))
## }

## ##
## ggcurve <- function(fun){
##     ggf = ggplot(data.frame(x=range(attr(fun,"support"))), aes(x)) +
##         stat_function(fun=fun, geom="line") +
##         tvcure_theme()
##     return(ggf)
## }

## ##
## plot.tvcureOld <- function(obj.tvcure,ci.level=.95,nrow=NULL,ncol=1,width=5,height=5){
##     obj = additive.tvcure(obj.tvcure)
##     z.alpha = qnorm(1-.5*(1-ci.level))
##     ##
##     ## depth.NULL <- function(x, ...) { browser(); 1; }
##     dev.new(width=10,height=5)
##     ## par(mfrow=c(1,2))
##     ggf0 = ggplot(data.frame(x=attr(obj$f0,"support")), aes(x)) +
##         stat_function(fun=obj$f0, geom="line") +
##         xlab("time") +
##         ylab(bquote(f[0](t))) +
##         tvcure_theme()
##     ggF0 = ggplot(data.frame(x=attr(obj$F0,"support")), aes(x)) +
##         stat_function(fun=obj$F0, geom="line") +
##         xlab("time") +
##         ylab(bquote(F[0](t))) +
##         tvcure_theme()
##     cat("ici-1\n")
##     ggf0
##     ggF0
##     grid.arrange(ggf0,ggF0,nrow=1,ncol=2,newpage=TRUE)
##     cat("ici-1b\n")
##     ## ##
##     ## addplot = function(f.list, fse.list, sub=1, ylim=NULL){
##     ##     J = length(f.list)
##     ##     ## Compute <ylim>
##     ##     if (is.null(ylim)){
##     ##         temp = c()
##     ##         for (j in 1:J) temp = c(temp,attr(f.list[[j]],"range"))
##     ##         ylim = range(temp)
##     ##     }
##     ##     if (is.null(nrow)) nrow = ceiling(J/ncol)
##     ##     plt = list()
##     ##     cnt = 0
##     ##     for (j in 1:J){ ## Loop over functions
##     ##         cnt = cnt + 1
##     ##         f = f.list[[j]] ; f.se = fse.list[[j]]
##     ##         lab = attr(fun,"label")
##     ##         ggf1 = ggcurve(fun=f) +
##     ##             ylim(ylim) +
##     ##             xlab(lab) +
##     ##             ylab(bquote(f[.(sub)](.(lab)))) +
##     ##             geom_function(fun=function(x) return(f(x) + .1*sin(x)), linetype="dashed") +
##     ##             geom_function(fun=function(x) return(f(x) - .1*sin(x)), linetype="dashed")
##     ##             ## geom_function(fun=function(x) return(f(x) + z.alpha*f.se(x)), linetype="dashed") +
##     ##             ## geom_function(fun=function(x) return(f(x) - z.alpha*f.se(x)), linetype="dashed")
##     ##             ## stat_function(fun=function(x) return(f(x) + z.alpha*f.se(x)), geom="line", linetype="dashed") +
##     ##             ## stat_function(fun=function(x) return(f(x) - z.alpha*f.se(x)), geom="line", linetype="dashed")
##     ##         plt[[cnt]] = ggf1
##     ##         if (cnt == (nrow*ncol) | j==J){ ## Plot as soon as the foreseen grid is full or the loop finished
##     ##             dev.new(width=width,height=height)
##     ##             if (j == J) nrow = ceiling(cnt/ncol)
##     ##             grid.arrange(grobs=plt,nrow=nrow,ncol=ncol) ##,newpage=TRUE)
##     ##             cnt = 0 ; plt = list()
##     ##         }
##     ##     }
##     ## }
##     ##
##     addplot <- function(x.grid, f.grid, fse.grid, sub=1, ED=NULL, ylim=NULL){
##         J = ncol(f.grid)
##         ## Pointwise CI
##         f.low = f.grid - z.alpha*fse.grid ## Pointwise credible interval lower limit
##         f.up  = f.grid + z.alpha*fse.grid ## Pointwise credible interval upper limit
##         ## Credible envelope
##         if (!is.null(ED)){
##             f.Low = f.Up = f.grid
##             for (j in 1:J){
##                 coef = sqrt(qchisq(ci.level, ED[j,1]))
##                 f.Low[,j] = f.grid[,j] - coef*fse.grid[,j] ## Credible envelope lower limit
##                 f.Up[,j]  = f.grid[,j] + coef*fse.grid[,j] ## Credible envelope upper limit
##             }
##         }
##         ## Compute <ylim>
##         if (is.null(ylim)){
##             if (is.null(ED)) ylim = range(cbind(f.low,f.up))
##             if (!is.null(ED)) ylim = range(cbind(f.low,f.up,f.Low,f.Up))
##         }
##         if (is.null(nrow)) nrow = ceiling(J/ncol)
##         plt = list()
##         cnt = 0
##         for (j in 1:J){ ## Loop over functions
##             cnt = cnt + 1
##             lab = colnames(f.grid)[j]
##             G = data.frame(x=x.grid[,j], y=f.grid[,j],
##                            ylow=f.low[,j], yup=f.up[,j])
##             if (!is.null(ED)){
##                 G$yLow = f.Low[,j] ; G$yUp = f.Up[,j]
##             }
##             col = "gray40"
##             ggf1 = ggplot(data = G, aes(x)) +
##                 geom_hline(yintercept = 0, color="grey", linetype="dotted") +
##                 geom_line(data = G, aes(x,y)) +
##                 geom_line(data = G, aes(x,ylow), linetype="dashed", color=col) +
##                 geom_line(data = G, aes(x,yup),  linetype="dashed", color=col) +
##                 ylim(ylim) +
##                 xlab(lab) +
##                 ylab(bquote(f[.(sub)](.(lab)))) +
##                 tvcure_theme()
##             if (!is.null(ED)){
##                 ggf1 = ggf1 +
##                 geom_line(data = G, aes(x,yLow), linetype="dotted", color=col) +
##                 geom_line(data = G, aes(x,yUp),  linetype="dotted", color=col)
##             }
##                 ## geom_function(fun=function(x) return(f(x) + z.alpha*f.se(x)), linetype="dashed") +
##                 ## geom_function(fun=function(x) return(f(x) - z.alpha*f.se(x)), linetype="dashed")
##                 ## stat_function(fun=function(x) return(f(x) + z.alpha*f.se(x)), geom="line", linetype="dashed") +
##                 ## stat_function(fun=function(x) return(f(x) - z.alpha*f.se(x)), geom="line", linetype="dashed")
##             plt[[cnt]] = ggf1
##             if (cnt == (nrow*ncol) | j==J){ ## Plot as soon as the foreseen grid is full or the loop finished
##                 dev.new(width=width,height=height)
##                 if (j == J) nrow = ceiling(cnt/ncol)
##                 cat("ici-2\n")
##                 grid.arrange(grobs=plt,nrow=nrow,ncol=ncol) ##,newpage=TRUE)
##                 cat("ici-2b\n")
##                 cnt = 0 ; plt = list()
##             }
##         }
##     } ## End addplot
##     ##
##     if (obj$J1 > 0){
##         addplot(x.grid=obj$x1.grid, f.grid=obj$f1.grid, fse.grid=obj$f1.grid.se, ED=obj$ED1, sub=1)
##     }
##     ##
##     if (obj$J2 > 0){
##         addplot(x.grid=obj$x2.grid, f.grid=obj$f2.grid, fse.grid=obj$f2.grid.se, ED=obj$ED2, sub=2)
##     }
##     ## if (obj$J1 > 0){
##     ##     addplot(f.list=obj$f1, fse.list=obj$f1.se, sub=1)
##     ## }
##     ## ##
##     ## if (obj$J2 > 0){
##     ##     addplot(f.list=obj$f2, fse.list=obj$f2.se, sub=2)
##     ## }
##     ## ##
##     ## if (obj$J1 > 0){
##     ##     ## dev.new(width=10,height=5)
##     ##     ylim = range(obj$f1.grid)
##     ##     if (is.null(nrow)) nrow = ceiling(J1/ncol)
##     ##     plt = list()
##     ##     cnt = 0
##     ##     J = obj$J1
##     ##     for (j in 1:J){
##     ##         cnt = cnt + 1
##     ##         fun = obj$f1[[j]]
##     ##         lab = attr(fun,"label")
##     ##         ## dev.new(width=5,height=5)
##     ##         ggf1 = ggcurve(fun) +
##     ##             ylim(ylim) +
##     ##             xlab(lab) +
##     ##             ylab(bquote(f[1](.(lab))))
##     ##         ## ggf1 = ggplot(data.frame(x=attr(fun,"support")), aes(x)) +
##     ##         ##     stat_function(fun=fun, geom="line") +
##     ##         ##     xlab(lab) +
##     ##         ##     ylab(bquote(f[1](.(lab)))) +
##     ##         ##     ylim(ylim) +
##     ##         ##     tvcure_theme()
##     ##         ## dev.new(width=5,height=5)
##     ##         ## curve(fun(x),xlim=attr(fun,"support"),xlab=lab,ylab=paste("f1.",lab,sep=""))
##     ##         ## plot(ggf1)
##     ##         ## plt [[j]] = ggf1
##     ##         plt [[cnt]] = ggf1
##     ##         if (cnt == (nrow*ncol) | j==J){ ## Plot as soon as the foreseen grid is full or the loop finished
##     ##             dev.new(width=width,height=height)
##     ##             if (j == J) nrow = ceiling(cnt/ncol)
##     ##             grid.arrange(grobs=plt,nrow=nrow,ncol=ncol) ##,newpage=TRUE)
##     ##             cnt = 0 ; plt = list()
##     ##         }
##     ##     }
##     ## plt [[3]] = ggf1
##     ## plt [[4]] = ggf1
##     ## grp = gl(ceiling(J1/(nrow*ncol),nrow*ncol,J1)
##     ## grid.arrange(grobs=plt,ncol=ncol, newpage=TRUE)
##     ## }
## }


## credibleRegionPlot = function(x, y, ylow, yup,
##                               ylim=NULL, xlab="x",ylab="y"){
##     if (is.null(ylim)) ylim = range(c(y,ylow,yup))
##     G = data.frame(x=x, y=y,
##                    ylow=ylow, yup=yup)
##     col = "gray40"
##     obj = ggplot(G, aes(x)) +
##                  geom_line(data = G, aes(x,y)) +
##                  geom_line(data = G, aes(x,ylow), linetype="dashed", color=col) +
##                  geom_line(data = G, aes(x,yup),  linetype="dashed", color=col) +
##                  ylim(ylim) +
##                  xlab(xlab) +
##                  ylab(ylab) +
##                  tvcure_theme()
##     return(obj)
## }

## ## Plot a single additive term
## addtermPlot <- function(j, x.grid, f.grid, fse.grid, sub=1, ED, ylim=NULL, ci.level=.95){
##     z.alpha = qnorm(1-.5*(1-ci.level))
##     J = ncol(f.grid)
##     ## Pointwise CI
##     f.low = f.grid - z.alpha*fse.grid ## Pointwise credible interval lower limit
##     f.up  = f.grid + z.alpha*fse.grid ## Pointwise credible interval upper limit
##     ## Credible envelope
##     if (!is.null(ED)){
##         f.Low = f.Up = f.grid
##         ## for (j in 1:J){
##             coef = sqrt(qchisq(ci.level, ED[j,1]))
##             f.Low[,j] = f.grid[,j] - coef*fse.grid[,j] ## Credible envelope lower limit
##             f.Up[,j]  = f.grid[,j] + coef*fse.grid[,j] ## Credible envelope upper limit
##         ## }
##     }
##     ## Compute <ylim>
##     if (is.null(ylim)){
##         if (is.null(ED)) ylim = range(cbind(f.low,f.up))
##         if (!is.null(ED)) ylim = range(cbind(f.low,f.up,f.Low,f.Up))
##     }
##     if (is.null(nrow)) nrow = ceiling(J/ncol)
##     plt = list()
##     cnt = 0
##     ## for (j in 1:J){ ## Loop over functions
##         cnt = cnt + 1
##         lab = colnames(f.grid)[j]
##         G = data.frame(x=x.grid[,j], y=f.grid[,j],
##                        ylow=f.low[,j], yup=f.up[,j])
##         if (!is.null(ED)){
##             G$yLow = f.Low[,j] ; G$yUp = f.Up[,j]
##         }
##         col = "gray40"
##         ggf1 = ggplot(data = G, aes(x)) +
##             geom_hline(yintercept = 0, color="grey", linetype="dotted") +
##             geom_line(data = G, aes(x,y)) +
##             geom_line(data = G, aes(x,ylow), linetype="dashed", color=col) +
##             geom_line(data = G, aes(x,yup),  linetype="dashed", color=col) +
##             ylim(ylim) +
##             xlab(lab) +
##             ylab(bquote(f[.(sub)](.(lab)))) +
##             tvcure_theme()
##         if (!is.null(ED)){
##             ggf1 = ggf1 +
##                 geom_line(data = G, aes(x,yLow), linetype="dotted", color=col) +
##                 geom_line(data = G, aes(x,yUp),  linetype="dotted", color=col)
##         }
##         ## plt[[cnt]] = ggf1
##         ## if (cnt == (nrow*ncol) | j==J){ ## Plot as soon as the foreseen grid is full or the loop finished
##         ##     ## dev.new(width=width,height=height)
##         ##     if (j == J) nrow = ceiling(cnt/ncol)
##         ##     grid.arrange(grobs=plt,nrow=nrow,ncol=ncol) ##,newpage=TRUE)
##         ##     cnt = 0 ; plt = list()
##         ## }
##     ## }
##     return(ggf1)
## } ## End addplot

## ## Get <xlim> and <ylim> of a ggplot object
## getxlim = function(obj) ggplot_build(obj)$layout$panel_params[[1]]$x.range
## getylim = function(obj) ggplot_build(obj)$layout$panel_params[[1]]$y.range
