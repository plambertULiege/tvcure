% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.tvcure.R
\name{predict.tvcure}
\alias{predict.tvcure}
\title{Predict method for tvcure model fits}
\usage{
\method{predict}{tvcure}(object, newdata, ci.level=.95, ...)
}
\arguments{
\item{object}{A \code{\link{tvcure.object}}.}

\item{newdata}{A data frame in which to look for the 'id' (distinguishing the different units), 'time' and covariate values for which 'predictions' should be made. Time values for a given 'id' should be a series of consecutive integers starting with 1. If \code{newdata$id} does not exist, then predictions are assumed to concern a single unit with consecutive time values starting with 1.}

\item{ci.level}{Credible level for the reported estimates. (Default: 0.95).}

\item{...}{additional generic arguments.}
}
\value{
A data frame containing, in addition to the optional \code{newdata} entries, the following elements:
\itemize{
\item \code{Hp} : Matrix containing estimates of the cumulative population hazard \eqn{H_p(t|x_{1:t})} with its credible interval bounds at time \eqn{t} given the history of covariates.
\item \code{lHp} : Matrix containing estimates of the log cumulative population hazard \eqn{\log H_p(t|x_{1:t})} with its standard error and credible interval bounds at time \eqn{t} given the history of covariates.
\item \code{se.lHp} : Vector containing the standard errors of the estimated log cumulative population hazard at time \eqn{t} given the history of covariates.
\item \code{hp} : Matrix containing estimates of the population hazard \eqn{h_p(t|x_{1:t})} with its credible interval bounds at time \eqn{t} given the history of covariates.
\item \code{lhp} : Matrix containing estimates of the log population hazard \eqn{\log h_p(t|x_{1:t})} with its standard error and credible interval bounds at time \eqn{t} given the history of covariates.
\item \code{se.lhp} : Vector containing the standard errors of the estimated log population hazard at time \eqn{t} given the history of covariates.
\item \code{Sp} : Matrix containing estimates of the population survival fuction \eqn{S_p(t|x_{1:t})=\exp(-H_p(t|x_{1:t}))} with its credible interval bounds at time \eqn{t} given the history of covariates.
\item \code{pcure} : Matrix containing estimates of the conditional cure probability of a unit still at tisk at time \eqn{t}, \eqn{P(T=+\infty|T>t,x=x_t)}, with its credible interval bounds at time \eqn{t} if covariates remain constant from time \eqn{t}.
\item \code{llpcure} : Matrix containing estimates of the conditional log-log cure probability of a unit still at tisk at time \eqn{t}, \eqn{\log(-\log P(T=+\infty|T>t,x=x_t))}, with its standard error and credible interval bounds at time \eqn{t} if covariates remain constant from time \eqn{t}.
\item \code{se.llpcure} : Vector containing the standard errors of the estimated conditional log-log cure probability of a unit still at tisk at time \eqn{t}, \eqn{\log(-\log P(T=+\infty|T>t,x=x_t))}, if covariates remain constant from time \eqn{t}.}
}
\description{
Predicted values based on a tvcure object.
}
\examples{
\donttest{
require(tvcure)
## Simulated data generation
beta = c(beta0=.4, beta1=-.2, beta2=.15) ; gam = c(gam1=.2, gam2=.2)
data = simulateTVcureData(n=500, seed=123, beta=beta, gam=gam,
                          RC.dist="exponential",mu.cens=550)$rawdata
## TVcure model fitting
tau.0 = 2.7 ; lambda1.0 = c(40,15) ; lambda2.0 = c(25,70) ## Optional
model = tvcure(~z1+z2+s(x1)+s(x2), ~z3+z4+s(x3)+s(x4), data=data,
               tau.0=tau.0, lambda1.0=lambda1.0, lambda2.0=lambda2.0)

## Covariate profiles for which 'predicted' values are requested
newdata = subset(data, id==1 | id==4)[,-3] ## Focus on units 1 & 4
pred = predict(model,newdata)

## Visualize the estimated population survival fns for units 1 & 4
## par(mfrow=c(1,2))
with(subset(pred,id==1), plotRegion(time,Sp,main="Id=1",
                              ylim=c(0,1),xlab="t",ylab="Sp(t)"))
with(subset(pred,id==4), plotRegion(time,Sp,main="Id=4",
                              ylim=c(0,1),xlab="t",ylab="Sp(t)"))
}
}
\references{
Lambert, P. and Kreyenfeld, M. (2025).
Time-varying exogenous covariates with frequently changing values in double additive cure survival model: an application to fertility.
\emph{Journal of the Royal Statistical Society, Series A}. <doi:10.1093/jrsssa/qnaf035>
}
\author{
Philippe Lambert \email{p.lambert@uliege.be}
}
