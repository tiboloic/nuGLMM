##
# LT 28/10/2018
#
# alpha nu-glmm, based on Bolker's glmmTMB
#

require(glmmTMB)

##' Fits nuglmm models to multiavriate abundance data
##'
##' @name nuglmm
##' @title Generalized LMM for multivariate abundance data
##' @param formula model specification
##' @param data dataframe
##' @param \dots optional additional arguments. Currently none are used in any
##' methods.
##' @return a fitted nuglmm model.
##' @keywords models
##' @examples
##' data(Tasmania, package = "mvabund")
##' m1 = nuglmm(abund ~ treatment * block, data=Tasmania, family="poisson")
##' m2 = nuglmm(abund ~ treatment * block, data=Tasmania, family="nbinom2")
##' anova(m1, m2)
##' @importFrom glmmTMB glmmTMB
##' @importFrom lme4 subbars findbars mkReTrms nobars
##' @importFrom methods is
##' @importFrom stats var getCall pchisq anova
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
  reTerms = lme4:::findbars(formula)
  
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

