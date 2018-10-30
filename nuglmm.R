##
# LT 28/10/2018
#
# alpha nu-glmm, based on Bolker's glmmTMB
#

#' Title
#'
#' @param formula 
#' @param data 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
nuglmm <- function(formula, data = NULL, ...) {
  
  # check LHS
  environment(formula) <- parent.frame()
  
  lhs = eval(formula[[2]], data, enclos=environment(formula))
  
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
  p = ncol(lhs)
  n = nrow(lhs)
  
  ## substitute LHS by first column of y
  f = formula(Y[,1]~1)
  # replace Y by LHS of formula
  f[[2]][[2]] = formula[[2]]
  # replace 1 by RHS of formula
  f[[3]] = formula[[3]]
  
  # deal with offset and weights
  
  # try a glmmTMB fit
  m = tryCatch(glmmTMB(f, data, doFit=FALSE, ...))
  
  .valid.family = c("gaussian", "poisson", "nbinom2", "binomial")
  
  if (!m$family$family %in% .valid.family)
    stop(paste0("This implementation of nuglmm accepts only the following families: ", .valid.family))
  
  fr = m$fr
  
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
  newTerms = list(formula)
  
  
  # treat existing random effects
  reTerms = lme4:::findbars(formula)
  
  for (ireTerm in seq_along(reTerms)) {
    reTerm = reTerms[[ireTerm]]
    newTerms[[ireTerm + 1]] = makeOp(makeOp(reTerm[[2]], makeOp(reTerm[[3]], quote(.spp), quote(`:`)), reTerm[[1]]), quote(`(`))
  }
  
  # build random effect term on fixed parameters across species
  newTerms[[length(newTerms) + 1]] = makeOp(makeOp(lme4:::nobars(formula)[[3]], quote(.spp), quote(`|`)), op=quote(`diag`))
  
  nuformula = Reduce(function(f,term) {glmmTMB:::RHSForm(f) <- makeOp(glmmTMB:::RHSForm(f), term, quote(`+`)); return(f)},newTerms) 
  
  ## try to call glmmTMB
  # don't fit, just get the structure
  m.struct = glmmTMB(nuformula, data=fr., doFit = FALSE, ...)
  
  # needed to get mu_predict in report
  m.struct$data.tmb$whichPredict = 1:nrow(fr.)
  
  # do the fit
  m = glmmTMB:::fitTMB(m.struct)
  
  m$n=nrow(fr)
  m$p = p
  m$TMBstruct = m.struct
  
  class(m) <- c("nuglmm", class(m))
  m
  
}

nuglmm.anova = function(m1, m2, method="", nboot = 100, ncpus = 1) {
  # LT 23/10/2018
  # first version of nuglmm anova using PITtrap !
  # m1 is the null model
  # m2 is the alternate
  
  # simulate random effects
  # simulating random effects in R is painfull if not resampling
  # just resample rows in R for alpha version ?
  
  # idea 1:
  # ensure I have the predicted mu in report, change the random effects and use report
  # example:
  ## get best fit parameters
  # ps = m2$obj$env$last.par.best
  ## set random params to resamp value
  # ps[random] = resampled Z * b
  ## get mu
  # mu = obj$report(ps)$mu_predict
  
  # FIXME control arguments
  # models must have class nuglmm
  # models must have same n, p, etc.
  # how to cheack for nestedness ?
  
  LR.obs = m1$fit$objective - m2$fit$objective
  
  # weaken convergence criteria
  m1$TMBstruct$control$optCtrl[["rel.tol"]] = 1e-3
  m2$TMBstruct$control$optCtrl[["rel.tol"]] = 1e-3
  m1$TMBstruct$control$optCtrl[["abs.tol"]] = 1e-3
  m2$TMBstruct$control$optCtrl[["abs.tol"]] = 1e-3
  
  # warm start from bestfit parameters
  m1$TMBstruct$parameters = m1$obj$env$parList()
  m2$TMBstruct$parameters = m2$obj$env$parList()
  
  fn = function(orig, iboot, m1, m2) {
    y.star = simulateone(m1, iboot)
    return(refitone(m1, m2, y.star))
  }
  
  options(warn = -1)
  sims = boot(1:m1$n, fn, nboot-1, parallel="snow", ncpus=ncpus, m1=m1, m2=m2)
  options(warn = 0)
  
  LR.star = c(LR.obs, sims$t)
  return(mean(LR.obs<=LR.star))
}

simulateone = function(m, iboot) {
  # iboot is a resampling of the rows
  
  # resample ranefs
  #Zu = m$obj$env$data$Z %*% m$obj$env$parList()$b
  #Zu.star = as.vector(matrix(Zu, nrow=n)[iboot,])
  #pars = m$obj$env$last.par.best
  # resampling is too complicate in case there is a random effect on a quantitative predictor
  # just do it conditional to the random effects or simulate the random effects
  # to simulate, could just put beta to 0 and simulate and get mu ?
  
  # response
  y = m$obj$env$data$yobs
  n = m$n
  
  # get mus: without ranefs from report, with ranefs, simulate  
  # simulate ranefs (parametric bootstrap)
  mus = m$obj$env$simulate()$mu_predict
  
  # resampling rows function
  resamp = function(pits) {
    as.vector(matrix(pits, nrow=n)[iboot,])
  }
  
  # get dispersion parameters if any
  # FIXME, should check if any dispersion info rather than look at names
  #if (any(names(m$fit$par)=="betad")) {
  phi = with(m$obj$env, exp(data$Xd %*% parList()$betad + ifelse(is.null(data$doffset), 0, data$doffset)))
  #}
  y.star = switch(m$modelInfo$family$family, 
                  "gaussian" = {qnorm(resamp(pnorm(y, mus, phi)), mus, phi)},
                  "poisson" = {qpois(resamp(pmin(runif(length(y), min = ppois(y-1, mus), max = ppois(y, mus)), 1-1e-8)), mus)},
                  "nbinom1" = {0},
                  "nbinom2" = {qnbinom(resamp(pmin(runif(length(y),
                                                         min = pnbinom(y-1, mu = mus, size = phi),
                                                         max = pnbinom(y, mu = mus, size = phi)), 1-1e-8)), mu = mus, size=phi)},
                  "binomial" = {rbinom(resamp(runif(length(y), min = pbinom(y-1, 1, mus), max = pbinom(y, 1, mus))), 1, mus)})
  
  
  y.star       
  
}

# refit
refitone = function(m1, m2, y.star) {
  
  # change data in TMBstruct 
  m1$TMBstruct$data.tmb$yobs = y.star
  m2$TMBstruct$data.tmb$yobs = y.star
  
  # weaken convergence criteria
  #m1$TMBstruct$control$optCtrl[["rel.tol"]] = 1e-3
  #m2$TMBstruct$control$optCtrl[["rel.tol"]] = 1e-3
  
  # warm start from bestfit parameters
  #m1$TMBstruct$parameters = m1$obj$env$parList()
  #m2$TMBstruct$parameters = m2$obj$env$parList()
  
  # possibly change control parameters to speed thinks up
  # refit
  m1.star = glmmTMB:::fitTMB(m1$TMBstruct)
  m2.star = glmmTMB:::fitTMB(m2$TMBstruct)
  
  # return log of likelihood ratio
  return(m1.star$fit$objective - m2.star$fit$objective)
}

