##
# LT 15/10/2018
#
# a revived version of nu-glmm, based on Bolker's glmmTMB
#


nuglmm <- function(formula, data = NULL, family = gaussian(), ziformula = ~0, dispformula = ~1,
                   weights = NULL, offset = NULL, contrasts = NULL, na.action = na.fail,
                   se=TRUE, verbose = FALSE, doFit = TRUE)
{
  ## edited copy-paste from glFormula
  
  ## edited copy-paste from glFormula
  ## glFormula <- function(formula, data=NULL, family = gaussian,
  ##                       subset, weights, na.action, offset,
  ##                       contrasts = NULL, mustart, etastart,
  ##                       control = glmerControl(), ...) {
  call <- mf <- mc <- match.call()
  
  if (is.character(family)) {
    if (family=="beta") {
      family <- "beta_family"
      warning("please use ",sQuote("beta_family()")," rather than ",
              sQuote("\"beta\"")," to specify a Beta-distributed response")
    }
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    ## call family with no arguments
    family <- family()
  }
  ## FIXME: what is this doing? call to a function that's not really
  ##  a family creation function?
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  fnames <- names(family)
  if (!all(c("family","link") %in% fnames))
    stop("'family' must contain at least 'family' and 'link' components")
  if (length(miss_comp <- setdiff(c("linkfun","variance"),fnames))>0) {
    warning("some components missing from ",sQuote("family"),
            ": downstream methods may fail")
  }
  if (grepl("^quasi", family$family))
    stop('"quasi" families cannot be used in glmmTMB')
  
  ## extract family and link information from family object
  link <- family$link
  
  ## lme4 function for warning about unused arguments in ...
  ## ignoreArgs <- c("start","verbose","devFunOnly",
  ##   "optimizer", "control", "nAGQ")
  ## l... <- list(...)
  ## l... <- l...[!names(l...) %in% ignoreArgs]
  ## do.call(checkArgs, c(list("glmer"), l...))
  
  # substitute evaluated versions
  ## FIXME: denv leftover from lme4, not defined yet
  
  environment(formula) <- parent.frame()
  call$formula <- mc$formula <- formula
  ## add offset-specified-as-argument to formula as + offset(...)
  ## need evaluate offset within envi
  if (!is.null(eval(substitute(offset),data,
                    enclos=environment(formula)))) {
    formula <- addForm0(formula,makeOp(substitute(offset),op=quote(offset)))
  }
  
  
  environment(ziformula) <- environment(formula)
  call$ziformula <- ziformula
  
  environment(dispformula) <- environment(formula)
  call$dispformula <- dispformula
  
  ## now work on evaluating model frame
  m <- match(c("data", "subset", "weights", "na.action", "offset"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  
  ## replace . in ziformula with conditional formula, ignoring offset
  if (glmmTMB:::inForm(ziformula,quote(.))) {
    ziformula <-
      update(glmmTMB:::RHSForm(glmmTMB:::drop.special2(formula),as.form=TRUE),
             ziformula)
  }
  
  ## want the model frame to contain the union of all variables
  ## used in any of the terms
  ## combine all formulas
  formList <- list(formula, ziformula, dispformula)
  for (i in seq_along(formList)) {
    f <- formList[[i]] ## abbreviate
    ## substitute "|" by "+"; drop specials
    f <- glmmTMB:::noSpecials(lme4:::subbars(f),delete=FALSE)
    formList[[i]] <- f
  }
  combForm <- do.call(glmmTMB:::addForm,formList)
  environment(combForm) <- environment(formula)
  ## model.frame.default looks for these objects in the environment
  ## of the *formula* (see 'extras', which is anything passed in ...),
  ## so they have to be put there ...
  for (i in c("weights", "offset")) {
    if (!eval(bquote(missing(x=.(i)))))
      assign(i, get(i, parent.frame()), environment(combForm))
  }
  
  mf$formula <- combForm
  fr <- eval(mf,envir=environment(formula),enclos=parent.frame())
  
  ## FIXME: throw an error *or* convert character to factor
  ## convert character vectors to factor (defensive)
  ## fr <- factorize(fr.form, fr, char.only = TRUE)
  ## store full, original formula & offset
  ## attr(fr,"formula") <- combForm  ## unnecessary?
  nobs <- nrow(fr)
  weights <- as.vector(model.weights(fr))
  
  #if(!is.null(weights) & !okWeights(family$family)) {
  #  stop("'weights' are not available for this family.")
  #}
  
  #if (is.null(weights)) weights <- rep(1,nobs)
  
  ## sanity checks (skipped!)
  ## wmsgNlev <- checkNlevels(reTrms$ flist, n=n, control, allow.n=TRUE)
  ## wmsgZdims <- checkZdims(reTrms$Ztlist, n=n, control, allow.n=TRUE)
  ## wmsgZrank <- checkZrank(reTrms$Zt, n=n, control, nonSmall=1e6, allow.n=TRUE)
  
  ## store info on location of response variable
  respCol <- attr(terms(fr), "response")
  names(respCol) <- names(fr)[respCol]
  
  ## extract response variable
  ## (name *must* be 'y' to match guts of family()$initialize
  y <- fr[,respCol]
  if (!is.matrix(y)) {
    stop("nuglmm needs a matrix-valued response")
  } else {
    if (ncol(y) < 5)
      stop("A minimumn of 5 columns (species) is required")
  }
  
  p = ncol(y)
  
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
  # don't fit, just get the TMBstruct
  m.struct = glmmTMB(nuformula, data=fr., family = family, doFit = FALSE)
  
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

