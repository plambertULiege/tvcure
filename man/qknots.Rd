% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qknots.R
\name{qknots}
\alias{qknots}
\title{Specification of the knots of a cubic B-spline basis for given data.}
\usage{
qknots(x, xmin=NULL, xmax=NULL,
       equid.knots = TRUE, pen.order=2, K=25)
}
\arguments{
\item{x}{data that should be supported by the knots of the B-spline basis.}

\item{xmin}{(Optional) minimum value for the knots.}

\item{xmax}{(Optional) maximum value for the knots.}

\item{equid.knots}{Logical indicating if equidistant knots are desired (Default: TRUE).}

\item{pen.order}{penalty order (if equid.knots = TRUE) (Default: 2).}

\item{K}{number of B-splines in the basis (Default: 25).}
}
\value{
a list containing the following elements:
\itemize{
\item \code{xmin} : minimum value of the knots.
\item \code{xmax} : maximum value of the knots.
\item \code{knots} : vector containing the knots: equidistant if \code{equid.knots} is TRUE, based on quantiles of \code{x} otherwise.
\item \code{Pd} : penalty matrix for the B-spline coefficients.
\item \code{pen.order} : penalty order for the P-spline model.
}
}
\description{
Specification of the knots of a cubic B-spline basis for given data.
}
\examples{
x = rnorm(100)
qknots(x)
}
\references{
Lambert, P. and Kreyenfeld, M. (2025).
Time-varying exogenous covariates with frequently changing values in double additive cure survival model: an application to fertility.
\emph{Journal of the Royal Statistical Society, Series A}. <doi:10.1093/jrsssa/qnaf035>
}
\author{
Philippe Lambert \email{p.lambert@uliege.be}
}
