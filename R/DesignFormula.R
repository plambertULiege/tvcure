## ---------------------------------------------------------------------------------------------------
## Philippe LAMBERT (ULiege, Oct 2018 ; updated to handle s(x,ref=val) in Oct 2022)(ULiege, Oct 2023)
## Email:  p.lambert@uliege.be
## Web: http://www.statsoc.ulg.ac.be
## ---------------------------------------------------------------------------------------------------
#' Internal function extracting design matrices from formulas in the tvcure function and computing penalty related matrices
#' @description Internal function extracting design matrices from formulas in the tvcure function and computing penalty related matrices.
#'
#' @param formula A formula describing the fixed effects and the additive terms in a regression model.
#' @param data A dataframe containing the data.
#' @param K Number of B-splines to describe an additive term.
#' @param pen.order Desired penalty order for the spline parameters in the additive terms.
#' @param knots.x (Optional) list of length J with the knots associated to each of the J additive terms. Automatically specified from the data by default.
#' @param n Number of units (Default: number of rows in the design matrix constructed from the formula and the data frame).
#' @param nointercept Logical indicating if the intercept should be set to zero (Default: FALSE).
#'
#' @return A list with
#' \itemize{
#' \item \code{J} : number of additive terms.
#' \item \code{K} : number of B-splines in a basis used to estimate an additive term.
#' \item \code{Z} : (n x nfixed) design matrix with fixed effects (including a first column of 1 if nointercept is FALSE).
#' \item \code{X} : (n x J) design matrix with the covariates involved in the additive terms.
#' \item \code{Xcal} : Z column-stacked with the J centered B-spline bases to yield the full design matrix (with column labels).
#' \item \code{nfixed} : number of fixed effect regression parameters.}
#'
#' If additive terms are specified in the formula, the following elements also appear:
#' \itemize{
#' \item \code{Bcal} : column-stacked matrix with the J centered B-spline bases.
#' \item \code{Bx} : list with J objects (one per additive term) including (B,Dd,Pd,K,cm).
#' \item \code{Pd.x, Dd.x} : penalty and difference penalty matrices applied on the spline parameters of an additive term.
#' \item \code{knots.x} : list of length J with the knots associated to each of the J additive terms.
#' \item \code{pen.order} : penalty order for the spline parameters in the additive terms.
#' \item \code{additive.lab} : labels for the columns in <Bcal> associated to the additive terms.
#' \item \code{lambda.lab} : labels for the penalty parameters.
#' \item \code{has.ref} : vector of J logicals indicating whether reference values were specified for a given additive term.
#' \item \code{ref.values} : specified reference values for the the J additive terms.
#' \item \code{cm.values} : list of length J with the values of the B-spline basis at the reference values specified (if any) for each of the additive terms.}
#'
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#' @references Lambert, P. and Kreyenfeld, M. (2024). Exogenous time-varying covariates in double additive cure survival model
#' with application to fertility. \emph{Journal of the Royal Statistical Society, Series A}, under review.
#'
#' @keywords internal
#'
DesignFormula = function(formula, data, K=10, pen.order=2, knots.x=NULL, n=NULL, nointercept=FALSE){
  if (!inherits(formula, "formula"))
    stop("Incorrect model formula")
  if ((formula=="~1")&(missing(data))){
    if (is.null(n)){
      cat("Model with only the intercept: the sample size <n> or a data frame should be provided !\n")
      return(NULL)
    }
    XX = model.matrix(~ 1, data = data.frame(rep(1,n)))
  } else {
    ## Extract design matrix
    if (missing(data)) {
      mf <- stats::model.frame(formula)
      XX <- stats::model.matrix(mf)
      if (nointercept){
          if (colnames(XX)[1] == "(Intercept)") XX = XX[,-1]
      }
    }
    else {
      mf <- stats::model.frame(formula, data = data)
      XX <- stats::model.matrix(mf, data = data)
      if (nointercept){
          if (colnames(XX)[1] == "(Intercept)") XX = XX[,-1,drop="FALSE"]
      }
    }
  }
  ## Identify additive terms
  smterms <- grepl("s(", colnames(XX), fixed = TRUE)
  X = subset(XX, select=smterms) ## Covariates with additive effects in the design matrix
  ## Reorder the design matrix to start with 'fixed effect' covariates
  idx = c(which(!smterms),which(smterms))
  XX <- subset(XX,select=idx) ## Reordering
  if (any(is.infinite(XX)))
    stop("Covariates contain Inf, NA or NaN values")
  J <- sum(smterms) ## Nbr of additive terms
  nfixed <- ncol(XX) - J
  n <- nrow(XX) ## Nbr of units
  ##
  if (nfixed == 0) Z = NULL
  else Z <- subset(XX, select=1:nfixed) ## 'Fixed effects' part of the design matrix
  ## if (ncol(Z) == 1)
  ##   colnames(Z) <- "(Intercept)"
  ## Some labels
  fixed.lab = colnames(Z) ## Labels of the fixed effects
  additive.lab = lambda.lab = Bcal = NULL
  ## Additives terms
  Bx = NULL ## B-spline for the additive terms
  if (J > 0){
      has.ref = rep(FALSE,J) ## Indicate if reference value for jth additive term
      ref.values = rep(NA,J) ## Possible covariate reference for jth additive term
      cm.values = vector(mode="list", length=J) ## B-splines values at reference for jth additive term
      K = floor(K) ## Make sure that this is an integer
      ## Number of knots
      nknots = K-1 ## Number of knots for the cubic B-splines in the basis
      if (!is.vector(nknots, mode = "numeric") || length(nknots) > 1 || is.na(nknots))
          stop("nknots must be numeric of length 1")
      if (K < 5 || K > 60)
          stop("K must be between 5 and 60")
      ## Penalty order
      pen.order <- floor(pen.order)
      if (pen.order < 1 || pen.order > 4)
          stop("Penalty order must be between 1 and 4")
      ## knots.x = NULL
      set.knots = is.null(knots.x)
      for (j in 1:J){  ## Loop over functional components in location
          xj = XX[,nfixed+j] ## Covariate for the jth additive term
          if (set.knots) knots.x[[j]] = seq(min(xj),max(xj),length=nknots) ## Knots
          colname = colnames(XX)[nfixed+j] ## Model specification for jth additive term in formula
          has.ref[j] = grepl("ref =", colname, fixed=TRUE) ## Detect specification of a reference value for the jth additive term
          if (has.ref[j]){
              txt = sub(")",'', sub(".*, ref = ",'',colname))
              if (txt == "NULL"){
                  ref.values[j] = NA
                  has.ref[j] = FALSE
              } else {
                  ref.values[j] = eval(parse(text=txt))
                  ## ref.values[j] = as.numeric(sub(")",'', sub(".*, ref = ",'',colname)))
                  cm.values[[j]] = head(c(Bsplines(ref.values[j], knots=knots.x[[j]])), -1)
              }
              ## Remove
              colnames(XX)[nfixed+j] = paste(sub(",.*",'',colname),")",sep='')
                  ## sub(",.*",'',sub("s\\(",'',colname))
          }
          ## B-spline matrix for the jth additive term
          Bx[[j]] = centeredBasis.gen(xj,knots=knots.x[[j]],cm=cm.values[[j]],pen.order)
          colnames(Bx[[j]]$B) = paste(colnames(XX)[nfixed+j],".",1:K,sep="")
    }
    ## Global design matrix
    Bcal = NULL ## Global design matrix for the B-splines of the additive terms
    for (j in 1:J) Bcal = cbind(Bcal,Bx[[j]]$B)
    ## Extra labels
    # additive.lab = colnames(XX)[-(1:nfixed)]
    if (nfixed > 0) additive.lab = unname(sapply(colnames(XX)[-(1:nfixed)], function(x) substring(x,3,nchar(x)-1)))
    else additive.lab = unname(sapply(colnames(XX), function(x) substring(x,3,nchar(x)-1)))
    lambda.lab = paste("lambda.",additive.lab,sep="")
  }
  if (J > 0) {
    Xcal = cbind(Z,Bcal) ## Global design matrix
  } else Xcal = Z
  ##
  Pd.x = Dd.x = NULL
  if (J > 0){
    Pd.x = Bx[[1]]$Pd ## Same penalty matrix for all functional components
    Dd.x = Bx[[1]]$Dd ## Same difference matrix for all functional components
  }
  ##
  ans = list(J=J,K=K,Z=Z,X=X,Xcal=Xcal,
             nfixed=nfixed)
  if (J > 0){
    ans = c(ans,
            list(
                Bcal=Bcal,
                Bx=Bx,Pd.x=Pd.x,Dd.x=Dd.x,knots.x=knots.x,
                pen.order=pen.order,additive.lab=additive.lab,lambda.lab=lambda.lab,
                has.ref=has.ref, ref.values=ref.values, cm.values=cm.values))
  }
  return(ans)
}
