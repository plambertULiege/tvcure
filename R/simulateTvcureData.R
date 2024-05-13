#' Simulation of survival data with a cure fraction and time-varying covariates.
#'
#' @description
#' Simulation of cure survival data in a counting process format with time-varying covariates.
#' The population hazard at time t  underlying the tvcure model is, for given covariate values,
#' \deqn{h_p(t|{\bf v}(t),\tilde{\bf v}(t)) = \mathrm{e}^{\eta_\vartheta({\bf v}(t))+\eta_F(\tilde{\bf v}(t))}
#' f_0(t)S_0(t)^{\exp(\eta_F(\tilde{\bf v}(t)))-1}}
#' with linear predictors
#' \deqn{\eta_\vartheta({\bf v}(t)) = \beta_0 + \beta_1 z_1(t) + \beta_2 z_2 + f_1(x_1(t)) + f_2(x_2(t))}
#' \deqn{\eta_F(\tilde{{\bf v}}(t)) = \gamma_1 z_3(t) +  \gamma_2 z_4 +  \tilde{f}_1(x_3(t)) + \tilde{f}_2(x_4(t))}
#' where \eqn{{\bf v}(t)=(z_1(t),z_2,x_1(t),x_2(t))}, \eqn{\tilde{\bf v}(t)=(z_3(t),z_4,x_3(t),x_4(t))},
#' with time-varying covariates \eqn{x_1(t)}, \eqn{x_3(t)} assumed identical and shared by the 2 submodels when
#' \code{shared.cov} is TRUE.
#'
#' The density \eqn{f_0(t)} governing the reference cumulative hazard dynamic is,
#' by default, a Weibull with shape parameter 2.65 and scale parameter 100,
#' ensuring that all susceptible units will experience the monitored event by time Tmax=300.
#'
#' The functions defining the additive terms are
#'  \deqn{f_1(x_1)= -.63 + .57\arctan(4x_1) ~;~f_2(x_2)= -.5 \cos(2\pi x_2)}
#'  \deqn{\tilde{f}_1(\tilde{x}_3) = .15 - .5 \cos\big(\pi(\tilde{x}_3-.75)\big)~;~
#'  \tilde{f}_2(\tilde{x}_4) = .8\big(\tilde{x}_4-.5\big)}
#'
#' Covariates are generated as follows:
#' \itemize{
#'  \item{ } \eqn{z_1(t), z_3(t)} are recentered time-varying Bernoulli(0.5) on \eqn{(0,T_{max})} ;
#'  \item{ } \eqn{z_2, z_4 \sim N(0,1)} ;
#'  \item{ } \eqn{x_1(t), x_2(t), x_3(t), x_4(t)} follow random cubic polynomial trajectories on \eqn{(0,T_{max})}.
#' }
#' More details can be found in Lambert and Kreyenfeld (2024).
#'
#' @usage simulateTVcureData(n, seed, Tmax=300,
#'        f0F0 = list(f0=function(x) dweibull(x, 2.65, 100),
#'                    F0=function(x) pweibull(x, 2.65, 100)),
#'        beta, gam, shared.cov=TRUE,
#'        RC.dist=c("uniform","exponential","Tmax"),
#'        tRC.min = 120, mu.cens=40, get.details=TRUE)
#' @param n Number of units.
#' @param seed Seed (integer) for the random TVcure data generator.
#' @param Tmax Maximum follow-up time after which a unit is considered cured in the absence of a previous event. (Default: 300).
#' @param f0F0 List of length 2 providing the density \eqn{f_0(t)} and associated cdf \eqn{F_0(t)} governing the bounded hazard dynamic on (0,Tmax), with \eqn{F_0}(Tmax)=1.0. (Default: f0F0 = list(f0=function(x) dweibull(x, 2.65, 100), F0=function(x) pweibull(x, 2.65, 100))).
#' @param beta 3-vector with the regression coefficients <beta> in the long-term (cure) survival (or quantum) submodel.
#' @param gam 2-vector with the regression coefficients <gamma> in the short-term (cure) survival (or timing) submodel.
#' @param shared.cov Logical indicating whether shared covariates for both submodels are assumed, with then \eqn{x_1(t)=x_3(t)}. (Default: TRUE).
#' @param RC.dist Right-censoring distribution: "uniform" (Uniform on (\code{tRC.min},\code{Tmax})),"exponential" (with mean \code{mu.cens}) or "Tmax" (when right-censoring occurs at Tmax)
#' @param tRC.min Minimum right-censoring time value if the right-censoring time distribution is Uniform. (Default: 120).
#' @param mu.cens Mean of the right-censoring time distribution if it is Exponential. (Default: 40).
#' @param get.details Logical indicating if a detailed data frame \code{rawdata} including the sequence of time-varying covariate values should also be reported. (Default: TRUE).
#'
#' @return A list with following elements:
#' \itemize{
#' \item{\code{seeds} : \verb{ }}{Seeds used to generate the data for each of the n units.}
#' \item{\code{tRC.min} : \verb{ }}{Minimum right-censoring time value if the right-censoring time distribution is Uniform.}
#' \item{\code{RC.dist} : \verb{ }}{Right-censoring distribution ("Uniform", "Exponential" or "Tmax").}
#' \item{\code{cure.rate} : \verb{ }}{Underlying proportion of cured units (i.e. without an observed event by \code{Tmax} if the follow-up is not interrupted by that time due to right-censoring).}
#' \item{\code{RC.rate} : \verb{ }}{Observed right-censoring rate.}
#' \item{\code{rawdata} : \verb{ }}{Data frame containing the generated data in a counting process format with the detailed follow-up for each unit until the event or right-censoring occurs:}
#'   \itemize{
#'   \item{\code{id} : \verb{ }}{Unit identificator for each row.}
#'   \item{\code{time} : \verb{ }}{Discrete observation times, starting at 1 for a given unit, until the end of its follow-up. The number of rows associated to a given unit corresponds to the follow-up duration.}
#'   \item{\code{event} : \verb{ }}{Event indicator (1 if it occured, 0 otherwise) for given unit at a given time.}
#'   \item{\code{z1, z2, z3, z4, x1, x2, x3, x4} : \verb{ }}{Covariate values for a given unit at a given time.}
#'   }
#' \item{\code{data.summary} : \verb{ }}{Data frame with n rows containing summarized information on the generated data for each unit:}
#'   \itemize{
#'   \item{\code{id} : \verb{ }}{Unit identificator (the ith row corresponding to the ith unit).}
#'   \item{\code{t.obs} : \verb{ }}{Observed event or right-censoring time.}
#'   \item{\code{delta} : \verb{ }}{Event indicator (1 if it occured, 0 otherwise).}
#'   \item{\code{t.true} : \verb{ }}{True (possibly unobserved) event time (Inf for a cured unit).}
#'   \item{\code{t.cens} : \verb{ }}{True (possibly unobserved) right-censoring time.}
#'   \item{\code{cured} : \verb{ }}{True (possibly unobserved) cure status.}
#'   }
#' \item{\code{parameters} : \verb{ }}{List containing the defining elements of the tvcure model:}
#'   \itemize{
#'   \item{\code{beta} : \verb{ }}{The regression parameters in the long-term survival (or quantum) submodel.}
#'   \item{\code{gam} : \verb{ }}{The regression parameters in the short-term survival (or timing) submodel.}
#'   \item{\code{f.theta} : \verb{ }}{A list of length 2 containing the functions defining the additive terms in the long-term survival (or quantum) submodel.}
#'   \item{\code{f.gam} : \verb{ }}{A list of length 2 containing the functions defining the additive terms in the short-term survival (or timing) submodel.}
#'   \item{\code{f0} : \verb{ }}{Density function governing the dynamic of the reference cumulative hazard on (0,Tmax).}
#'   \item{\code{F0} : \verb{ }}{CDF governing the dynamic of the reference cumulative hazard on (0,Tmax).}
#'   }
#' \item{\code{call} : \verb{ }}{Function call.}
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
#' ## Regression parameters
#' beta = c(beta0=.4, beta1=-.2, beta2=.15) ##  beta0 tunes the cure rate
#' gam = c(gam1=.2, gam2=.2)
#' ## Data simulation
#' temp = simulateTVcureData(n=500, seed=123, beta=beta, gam=gam,
#'                           RC.dist="exponential",mu.cens=550)
#' head(temp$rawdata) ## Overview of the simulated raw data
#' head(temp$data.summary) ## Overview of the summarized data
#' with(temp, c(cure.rate=cure.rate,RC.rate=RC.rate)) ## Cure and RC rates
#'
simulateTVcureData = function(n, seed, Tmax=300,
                        f0F0 = list(f0=function(x) dweibull(x, 2.65, 100), F0=function(x) pweibull(x, 2.65, 100)),
                        beta, gam, shared.cov=TRUE, RC.dist=c("uniform","exponential","Tmax"),
                        tRC.min = 120, mu.cens=40, get.details=TRUE){
    cl <- match.call()
    ## --------------------
    ## Covariate generation
    ## --------------------
    ##
    ## Generate a random polynomial function on (t1,t2) with values in (y1,y2)
    ## ----------------------------------------------------------------------
    ## Quadratic
    poly2.gen = function(t1=0,t2=1,y1=0,y2=1){
        t = c(t1,mean(c(t1,t2)),t2)
        y = y1 + (y2-y1)*.9*runif(length(t))
        coef = solve(cbind(1, t, t^2), y)
        fun = function(t) pmax(y1,pmin(y2,c(cbind(1, t, t^2) %*% coef)))
        return(list(t=t,y=y,coef=coef,fun=fun))
    } ## End of poly2.gen
    ##
    ## Cubic
    poly3.gen = function(t1=0,t2=1,y1=0,y2=1){
        t = seq(t1,t2,length=4)
        y = y1 + (y2-y1)*.9*runif(length(t))
        coef = solve(cbind(1, t, t^2, t^3), y)
        fun = function(t) pmax(y1,pmin(y2,c(cbind(1, t, t^2, t^3) %*% coef)))
        return(list(t=t,y=y,coef=coef,fun=fun))
    } ## End of poly3.gen
    ##
    ## Quartic
    poly4.gen = function(t1=0,t2=1,y1=0,y2=1){
        t = seq(t1,t2,length=5)
        y = y1 + (y2-y1)*.9*runif(length(t))
        coef = solve(cbind(1, t, t^2, t^3, t^4), y)
        fun = function(t) pmax(y1,pmin(y2,c(cbind(1, t, t^2, t^3, t^4) %*% coef)))
        return(list(t=t,y=y,coef=coef,fun=fun))
    } ## End of poly4.gen
    ##
    ## Generate a binary step-function with <nvals> alternating values
    ## on a regular grid on (t1,t2)
    ## ---------------------------------------------------------------
    binary.gen = function(t1=0,t2=1,nvals,alternate=TRUE,regular=TRUE){
        eps = 1e-6
        cuts = seq(t1,t2,length=nvals+1) ## Regular changepoints
        if (!regular) cuts = c(t1,sort(runif(nvals-1,t1,t2)),t2+eps) ## Irregular Changepoints
        cuts[nvals+1] = t2 + eps
        if (alternate){
            vals = (0:(nvals-1)) %% 2 ## Sequence of 0-1 starting with 0
            if (runif(1)>=.5) vals <- c(1,head(vals,nvals-1)) ## Randomly start with 0 or 1
        } else {
            vals = rbinom(nvals,1,.5)
        }
        vals = vals-.5 ## Recentering of the binary covariate: (0,1) --> (-.5,.5)
        fun = function(x) return( vals[cut(x,cuts,include.lowest=TRUE)]) ## Step function
        return(list(cuts=cuts,vals=vals,fun=fun))
    }
    ## End of binary.gen
    ##
    ## Generate constant (z2,z4) and time-varying (z1t,z3t,x1t,x2t,x3t) covariates
    ## ---------------------------------------------------------------------------
    cov.gen = function(seed){
        if (!missing(seed)) set.seed(seed)
        tmin = 0 ; tmax = 300
        ans = list(
            seed=seed,
            tmin=tmin, tmax=tmax,
            ## Time-varying binary covariates
            z1t = binary.gen(t1=tmin,t2=tmax,nvals=1+rpois(1,5),alternate=FALSE,regular=TRUE)$fun,
            z3t = binary.gen(t1=tmin,t2=tmax,nvals=1+rpois(1,5),alternate=FALSE,regular=TRUE)$fun,
            ## Continuous constant covariates
            z2 = rnorm(1), z4 = rnorm(1),
            ## Continuous time-varying covariates
            x1t = poly3.gen(t1=tmin,t2=tmax,y1=0,y2=1.5)$fun,
            x2t = poly3.gen(t1=tmin,t2=tmax,y1=0,y2=1)$fun,
            x4t = poly3.gen(t1=tmin,t2=tmax,y1=0,y2=1)$fun
        )
        if (!shared.cov) ans$x3t = poly3.gen(t1=tmin,t2=tmax,y1=0,y2=1.5)$fun
        return(ans)
    }
    ## End of cov.gen
    ##
    ## ----------------------------------------------------------------
    ## Log-hazard for single unit based on covariate daily trajectories
    ## ----------------------------------------------------------------
    loghaz = function(t,Xt,beta,gam,...){
        ## Baseline (standardized) hazard
        ## ------------------------------
        ## Cure suvival model: Sp(t) = exp(-theta(x) F(t|x))
        ##  with theta(x) = exp(eta1(x))  &  1-F(t|x) = (1-F0(t))^exp(eta2(x))
        ##
        ## f0 = function(x) dweibull(x, 2.65, 100)
        ## F0 = function(x) pweibull(x, 2.65, 100)
        f0 = f0F0[[1]] ; F0 = f0F0[[2]]
        f0.t = f0(t) ; F0.t = F0(t)
        ## Covariates for cure prob ("long-term survival")
        ## -----------------------------------------------
        ## ... Linear effect covariates
        z1t = Xt$z1t(t) ## TV binary
        z2 = Xt$z2 ## rnorm(n,0,1)
        ## ... Additive
        x1t = Xt$x1t(t) ## poly3 on (0,1.5) ## Generate random (fixed) income
        x2t = Xt$x2t(t) ## poly3 on (0,1)
        ## Covariates for timing submodel ("short-term survival")
        ## ------------------------------------------------------
        ## ... Fixed cov
        z3t = Xt$z3t(t) ## TV binary
        z4 = Xt$z4 ## rnorm(n,0,1)
        ## ... Additive
        if (shared.cov){
            x3t = x1t ## Shared covariate (income)
        } else {
            x3t = Xt$x3t(t)
        }
        x4t = Xt$x4t(t) ## poly3 on (0,1)
        ##
        ## Regression parameters & Additive terms for cure prob
        ## ----------------------------------------------------
        beta0 = beta[1] ; beta1 = beta[2] ; beta2 = beta[3]
        f1.theta = function(x) -.63+.57*atan(4*x)
        f2.theta = function(x) -.5*cos(2*pi*x) ## Extra artifical additive term with x in (0,1)
        ## Regression parameters & Additive terms for timing submodel
        ## ----------------------------------------------------------
        gam1 = gam[1] ; gam2 = gam[2]
        f1.gam = function(x) -(.5*cos(pi*(x-.75)) -.15)
        f2.gam = function(x) .8*(x-.5)
        ##
        ## Log-hazard function
        ## -------------------
        eta.v = beta0 + beta1 * z1t + beta2 * z2 + f1.theta(x1t) + f2.theta(x2t)
        eta.F =          gam1 * z3t +  gam2 * z4 +   f1.gam(x3t) + f2.gam(x4t)
        ##
        loghaz = eta.v + eta.F + log(f0.t) + (exp(eta.F)-1) * log(1-F0.t)
        ##
        data = data.frame(t=t,
                        loghaz=loghaz,
                        z1=z1t, z2=z2,
                        x1=x1t, x2=x2t,
                        z3=z3t, z4=z4,
                        x3=x3t, x4=x4t)
        ans = list(data=data,
                   beta=beta,
                   gam=gam,
                   f.theta=c(f1=f1.theta,f2=f2.theta),
                   f.gam=c(f1=f1.gam,f2=f2.gam),
                   f0=f0, F0=F0)
        return(ans)
    } ## End of loghaz
    ##
    ## -------------------------------------------------------
    ## Data generation for single unit given covariate history
    ## -------------------------------------------------------
    simulSingleSurv = function(seed,
                               Xt,
                               beta,gam,
                               mu.cens = 40, tRC.min = 120, Tmax=300,
                               RC.dist=c("uniform","exponential","Tmax"),
                               get.details=FALSE){
        if (!missing(seed)) set.seed(seed)
        ## Time grid
        T = 300
        t.grid = 1:T ; dt = 1
        t.mid = t.grid - .5*dt
        ## Log-hazard on time grid for given covariate history
        obj.lht = loghaz(t.mid, Xt=Xt, beta=beta, gam=gam)
        lht = obj.lht$data$loghaz
        ## ht, St, Ft, Fmax
        ht = exp(lht)
        Ht = cumsum(ht*dt)
        St = exp(-Ht)
        Ft = 1-St
        Fmax = tail(Ft,1) ## Maximum cdf value
        ## Time-to-event data generation
        u = runif(1) ## Random number generator
        ## Event time (Inf if no event)
        t.true = ifelse(u <= Fmax, min(t.grid[Ft>=u]), Inf)
        cured = is.infinite(t.true)
        ##
        switch(RC.dist, ## Right-censoring time
               "uniform" = { ## Uniform
                   t.cens = sample(tRC.min:(T-1),1,replace=TRUE)
               },
               "exponential" = { ## Exponential
                   t.cens = pmin(ceiling(rexp(1,1/mu.cens)), T)
               },
               "Tmax" = { # no RC (in addition to cure and max follow-up to T)
                   t.cens = Tmax
               }
               )
        delta = 0 + (t.true < t.cens) ## Final event indicator
        t.obs = ifelse(delta==1, t.true, t.cens) ## Reported time
        idx = which(t.grid <= t.obs) ## Index of observed t.grid values
        event = 0*t.grid ; event[length(idx)] = delta
        ##
        if (get.details){ ## Include lht,ht,Ht,St,Ft
            data = data.frame(t.mid=t.mid,
                            time=t.grid,
                            event=event,
                            lht=lht, ht=ht, Ht=Ht,
                            St=St, Ft=Ft)[idx,]
        } else {
            data = data.frame(time=t.grid,
                            event=event)[idx,]
        }
        data = cbind(data, obj.lht$data[idx,-(1:2)]) ## Add covariate values
        ##
        ans = list(data=data,
                   cured=cured, T=T,
                   Fmax=Fmax,
                   u=u, t.true=t.true,
                   t.obs=t.obs, delta=delta,
                   t.cens=t.cens,
                   parameters = obj.lht[-1]
                   )
        return(ans)
    } ## End of simulSingleSurv
    ##
    ## ------------------------------------------------------------------------
    ## MAIN function: data generation for n units and their covariate histories
    ## ------------------------------------------------------------------------
    RC.dist = match.arg(RC.dist)
    ##
    set.seed(seed)
    seeds = sample(1:(2^19),n,replace=FALSE) ##seed + (1:n) - 1
    donnees = list()
    rawdata = c()
    t.true = t.obs = delta = t.cens= cured = numeric(n)
    for (i in 1:n){
        Xt = cov.gen(seed=seeds[i])
        obj.i = simulSingleSurv(seed=seeds[i],Xt=Xt,
                                beta=beta, gam=gam,
                                tRC.min=tRC.min, mu.cens=mu.cens, Tmax=Tmax, RC.dist=RC.dist)
        t.true[i] = obj.i$t.true
        t.obs[i]  = obj.i$t.obs
        delta[i]  = obj.i$delta
        t.cens[i] = obj.i$t.cens
        cured[i]  = obj.i$cured
        ## Data frame in a person-month format
        id = rep(i, obj.i$t.obs)
        data = cbind(id=id, obj.i$data) ## Data frame for unit i
        ## Data frame combining all units
        if (get.details) rawdata = rbind(rawdata,data)
    }
    ## Also return a data frame in a summarized format (1 line per unit)
    data.summary = data.frame(id=1:n,
                            t.obs=t.obs, delta=delta,
                            t.true=t.true, t.cens=t.cens,
                            cured=cured)
    ##
    if (!get.details){
        ans = list(data.summary=data.summary,
                   n=n,beta=beta,gam=gam,
                   cure.rate = sum(cured)/n,
                   RC.rate = sum(t.obs==t.cens)/n)
        return(ans)
    }
    ##
    donnees = list(seeds=seeds,
                   tRC.min=tRC.min, RC.dist=RC.dist,
                   cure.rate = sum(cured)/n,
                   RC.rate = sum(t.obs==t.cens)/n,
                   rawdata=rawdata,
                   data.summary=data.summary,
                   parameters=obj.i$parameters,
                   call=cl)
    return(donnees)
}
