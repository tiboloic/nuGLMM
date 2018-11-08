##
# LT 28/10/2018
#
# alpha nu-glmm, based on Bolker's glmmTMB
#

##' @title Generalized linear mixed effect models for multivariate abundance data
##' @description Fit glmms to multivariate abundance data using glmmTMB
##' @param formula combined fixed and random effects formula, following lme4 and glmmTMB
##'     syntax. The left hand side must evaluate to an abundance matrix with species in columns.
##' @param data data frame
##' @param family a family function, a character string naming a family function, or the result of a call to a family function family (variance/link function) information; see \code{\link{family}} for generic discussion of families or \code{\link{family_glmmTMB}} for details of \code{glmmTMB}-specific families.
##'     Only \code{gaussian()}, \code{binomial()}, \code{poisson()} and \code{nbinom2()} are compatible with nuglmm.
##' @param ziformula 
##' Zero-inflation is not yet available in nuglmm.
##' @param dispformula a \emph{one-sided} formula for dispersion containing only fixed effects: the
##'     default \code{~1} specifies the standard dispersion given any family.
##'     The argument is ignored for families that do not have a dispersion parameter.
##'     For an explanation of the dispersion parameter for each family, see (\code{\link{sigma}}).
##'     The dispersion model uses a log link. 
##'     In Gaussian mixed models, \code{dispformula=~0} fixes the parameter to be 0, forcing variance into the random effects.
##' @param weights weights, as in \code{glm}. Not automatically scaled to have sum 1.
##' @param offset offset for conditional model (only)
##' @param contrasts an optional list, e.g. \code{list(fac1="contr.sum")}. See the \code{contrasts.arg} of \code{\link{model.matrix.default}}.
##' @param na.action how to handle missing values (see \code{\link{na.action}} and \code{\link{model.frame}}); from \code{\link{lm}}, \dQuote{The default is set by the \code{\link{na.action}} setting of \code{\link{options}}, and is \code{\link{na.fail}} if that is unset.  The \sQuote{factory-fresh} default is \code{\link{na.omit}}.}
##' @return a fitted nuglmm model.
##' @details
##' \itemize{
##' \item binomial models with more than one trial (i.e., not binary/Bernoulli)
##' must be specified in the form \code{prob ~ ..., weights = N}.
##' }
##' @references
##' \itemize{
##' \item Warton, David & Thibaut, Loïc & Wang, Yi. (2017). The PIT-trap—A “model-free” bootstrap procedure for inference about regression models with discrete, multivariate responses. PLoS ONE. 12. 10.1371/journal.pone.0181790. 
##' \item Millar, Russell B. Maximum Likelihood Estimation and Inference: With Examples in R, SAS and ADMB. John Wiley & Sons, 2011.
##' }

##' @keywords models
##' @examples
##' data(Tasmania, package = "mvabund")
##' m1 = nuglmm(abund ~ treatment, data=Tasmania, family=nbinom2())
##' m2 = nuglmm(abund ~ treatment + (1|block), data=Tasmania, family=nbinom2())
##' anova(m1, m2)
##' @importFrom glmmTMB glmmTMB
##' @importFrom lme4 findbars
##' @importFrom methods is
##' @importFrom stats gaussian binomial poisson nlminb as.formula terms model.weights
##' @export
nuglmm <- function(formula, data = NULL, offset=NULL, weights = NULL, ...) {
  
  # save call
  cl <- match.call()
  
  # check LHS
  environment(formula) <- parent.frame()
  
  lhs <- eval(formula[[2]], data, enclos=environment(formula))
  
  ## need evaluate offset within envi
  if (!is.null(lhs)) {
    if (!is.matrix(lhs)) {
      stop("nuglmm needs a matrix-valued response")
    } else {
      if (ncol(lhs) < 5)
        stop("A minimumn of 5 columns (species) is required")
    }
  } else stop("Cannot evaluate response")
  
  # dimensions
  p <- ncol(lhs)
  n <- nrow(lhs)
  
  ## substitute LHS by first column of y
  f <- formula(Y[,1]~1)
  # replace Y by LHS of formula
  f[[2]][[2]] <- formula[[2]]
  # replace 1 by RHS of formula
  f[[3]] <- formula[[3]]
  
  
  # call glmmTMB to check arguments and get frame
  m <- tryCatch(glmmTMB(f, data, doFit=FALSE, ...), error = function(e) e)
  if ("error" %in% class(m))
    stop(m$message)
  
  .valid.family <- c("gaussian", "poisson", "nbinom2", "binomial")
  
  if (!m$family$family %in% .valid.family)
    stop(paste0("This implementation of nuglmm accepts only the following families: ", .valid.family))
  
  # retrieve frame
  fr <- m$fr
  
  respCol <- attr(terms(fr), "response")
  
  # restore original name for LHS
  names(respCol) <- deparse(formula[[2]])
  
  ## response variable
  ## (name *must* be 'y' to match guts of family()$initialize
  y <- lhs
  
  # modify model frame
  #  stack predictors and factors
  fr. = as.data.frame(lapply(fr[-respCol], rep, p))
  fr.[names(respCol)] = as.vector(y)
  fr.[".spp"]= gl(p, nrow(fr))
  
  ## start terms manipulation
  nuformula = formula
  
  # we need the new makeOp
  makeOp <- function(x,y,op=NULL) {
    if (is.null(op) || missing(y)) {  ## unary
      if (is.null(op)) {
        substitute(OP(X),list(X=x,OP=y))
      } else {
        substitute(OP(X),list(X=x,OP=op))
      }
    } else substitute(OP(X,Y), list(X=x,OP=op,Y=y))
  }
  
  # list of terms to be added, initialized with original formula
  newTerms <- list(formula)
  
  
  # treat existing random effects
  reTerms = lme4::findbars(formula)
  
  for (ireTerm in seq_along(reTerms)) {
    reTerm = reTerms[[ireTerm]]
    newTerms[[ireTerm + 1]] = makeOp(makeOp(reTerm[[2]], makeOp(reTerm[[3]], quote(.spp), quote(`:`)), reTerm[[1]]), quote(`(`))
  }
  
  # build random effect term on fixed parameters across species
  newTerms[[length(newTerms) + 1]] <- makeOp(makeOp(lme4:::nobars(formula)[[3]], quote(.spp), quote(`|`)), op=quote(`diag`))
  
  nuformula <- Reduce(function(f,term) {glmmTMB:::RHSForm(f) <- makeOp(glmmTMB:::RHSForm(f), term, quote(`+`)); return(f)},newTerms) 
  
  ## deal with offset and weights
  if (!missing(offset)) {
    if (is.matrix(offset)) {
      if (all(dim(offset) == c(n,p)))
        offset <- as.vector(offset)
      else
        stop("Incorrect dimension of offset matrix")
    }
  }
  
  if (!missing(weights)) {
    if (is.matrix(weights)) {
      if (all(dim(weights) == c(n,p)))
        weights <- as.vector(weights)
      else
        stop("Incorrect dimension of weight matrix")
    }
  }
  
  ## try to call glmmTMB
  # don't fit, just get the structure
  m.struct <- glmmTMB(nuformula, data=fr., offset = offset, weights = weights, doFit = FALSE, ...)
  
  # needed to get mu_predict in report
  m.struct$data.tmb$whichPredict <- 1:nrow(fr.)
  
  # do the fit
  m <- glmmTMB:::fitTMB(m.struct)
  
  m$n <- nrow(fr)
  m$p <- p
  m$TMBstruct <- m.struct
  
  class(m) <- c("nuglmm", class(m))
  m
  
}

