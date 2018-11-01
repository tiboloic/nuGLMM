##
# LT 31/10/2018
#
# methods for nuglmm objects
#

# copy pasted and edited from glmmTMB
# add parameter nboot,
# if nboot > 0, perform PITtrap

##' @importFrom methods is
##' @importFrom stats var getCall pchisq anova
##' @export
anova.nuglmm <- function (object, ..., nboot=0, model.names = NULL) 
{
  mCall <- match.call(expand.dots = TRUE)
  dots <- list(...)
  .sapply <- function(L, FUN, ...) unlist(lapply(L, FUN, ...))

  ## detect multiple models, i.e. models in ...
  modp <- as.logical(vapply(dots, is, NA, "nuglmm"))
  if (any(modp)) {
    mods <- c(list(object), dots[modp])
    nobs.vec <- vapply(mods, nobs, 1L)
    if (var(nobs.vec) > 0) 
      stop("models were not all fitted to the same size of dataset")
    if (is.null(mNms <- model.names)) 
      mNms <- vapply(as.list(mCall)[c(FALSE, TRUE, modp)], 
                     glmmTMB:::safeDeparse, "")
    if (any(duplicated(mNms))) {
      warning("failed to find unique model names, assigning generic names")
      mNms <- paste0("MODEL", seq_along(mNms))
    }
    if (length(mNms) != length(mods)) 
      stop("model names vector and model list have different lengths")
    names(mods) <- sub("@env$", "", mNms)

    llks <- lapply(mods, logLik)
    ii <- order(Df <- vapply(llks, attr, FUN.VALUE = numeric(1), 
                             "df"))
    mods <- mods[ii]
    llks <- llks[ii]
    Df <- Df[ii]
    calls <- lapply(mods, getCall)
    data <- lapply(calls, `[[`, "data")
    if (!all(vapply(data, identical, NA, data[[1]]))) 
      stop("all models must be fit to the same data object")
    header <- paste("Data:", glmmTMB:::abbrDeparse(data[[1]]))
    subset <- lapply(calls, `[[`, "subset")
    if (!all(vapply(subset, identical, NA, subset[[1]]))) 
      stop("all models must use the same subset")
    if (!is.null(subset[[1]])) 
      header <- c(header, paste("Subset:", glmmTMB:::abbrDeparse(subset[[1]])))
    llk <- unlist(llks)
    chisq <- 2 * pmax(0, c(NA, diff(llk)))
    dfChisq <- c(NA, diff(Df))

    if (nboot>0) {
      null.dists <- matrix(0, nboot+1, length(mods)-1)
      for (imod in 1:(length(mods)-1))
        null.dists[,imod] <- pittrap(mods[[imod]], mods[[imod + 1]], nboot, ...)
      val <- data.frame(Df = Df, AIC = .sapply(llks, AIC), 
                        BIC = .sapply(llks, BIC), logLik = llk, deviance = -2 * 
                          llk, LR = c(NA,null.dists[1,]),# `LRT Df` = dfChisq,
                        `Pr(>LR)` = c(NA, apply(null.dists, 2, function(null.dist) mean(null.dist[1]<=null.dist))), row.names = names(mods), 
                        check.names = FALSE)
      header <- c(paste0("Null distribution of likelihood ratio estimated by ", nboot, " PITtrap simulations"), header)
      
    } else {
    val <- data.frame(Df = Df, AIC = .sapply(llks, AIC), 
                      BIC = .sapply(llks, BIC), logLik = llk, deviance = -2 * 
                        llk, Chisq = chisq, `Chi Df` = dfChisq, `Pr(>Chisq)` = pchisq(chisq, 
                                                                                      dfChisq, lower.tail = FALSE), row.names = names(mods), 
                      check.names = FALSE)
    }
    class(val) <- c("anova", class(val))
    forms <- lapply(lapply(calls, `[[`, "formula"), deparse)
    ziforms <- lapply(lapply(calls, `[[`, "ziformula"), deparse)
    dispforms <- lapply(lapply(calls, `[[`, "dispformula"), deparse)
    #FIXME only output nontrivial ziforms and dispforms
    structure(val, heading = c(header, "Models:", 
                               paste(paste(paste(rep(names(mods), times = lengths(forms)), unlist(forms), sep = ": "),
                                           unlist(ziforms), sep=", zi="),
                                     unlist(dispforms), sep=", disp=")))
  } else stop("no single-model anova() method for nuglmm")
}

# build null distribution of likelihood ratio under PITtrap
# m1 is the null model, m2 is the alternative model, nboot is the number of simulations
pittrap <- function (m1, m2, nboot, ..., parallel = c("no", "multicore", "snow"),
                     ncpus = getOption("nuglmm.ncpus"),
                     cl = NULL) {
  
  # sort out parallel stuff
  # inspired by boot
  if (missing(ncpus))
    ncpus <- parallel::detectCores()
  if (missing(parallel)) {
    parallel <- getOption("nuglmm.parallel")
    if (is.null(parallel)) {
      # check number of cores
      parallel <- ifelse(parallel::detectCores()>1, "multicore", "no")
    }
  }
  parallel <- match.arg(parallel)
  have_mc <- have_snow <- FALSE
  if (parallel != "no" && ncpus > 1L) {
    if (parallel == "multicore") {
      have_mc <- .Platform$OS.type != "windows"
      have_snow <- .Platform$OS.type == "windows"
    }
    else if (parallel == "snow") 
      have_snow <- TRUE
    if (!have_mc && !have_snow) 
      ncpus <- 1L
    #loadNamespace("parallel")
  }
  
  # observed statistics
  LR.obs <- m1$fit$objective - m2$fit$objective
  
  # warm start from bestfit parameters
  m1$TMBstruct$parameters <- m1$obj$env$parList()
  m2$TMBstruct$parameters <- m2$obj$env$parList()
  
  # FIX ME: do parallel stuff
  
  if (have_snow)
    FUN = function(i, m1, m2) {
      # when nuglmm is released as a package, require(nuglmm) should be enough
      require(glmmTMB)
      source("Methods.R")
      y.star = simulateone(m1)
      return(refitone(m1, m2, y.star))
    }
  else
    FUN = function(i) {
      y.star = simulateone(m1)
      return(refitone(m1, m2, y.star))
    }
  
  if (ncpus > 1) {
    if (have_mc) {
      L <- parallel::mclapply(1:nboot, FUN, mc.cores = ncpus)
    } else if (have_snow) {
      if (is.null(cl)) {
        ## start cluster
        new_cl <- TRUE
        cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
      }
      ## run
      L <- parallel::clusterApply(cl, 1:nboot, FUN, m1=m1, m2=m2)
      if (new_cl) {
        ## stop cluster
        parallel::stopCluster(cl)
      }
    }} else { ## non-parallel
      L <- lapply(1:nboot, FUN)
    }
  
  L <- unlist(L)

  #options(warn = -1)
  #sims = boot(1:m1$n, fn, nboot, parallel="snow", ncpus=2, m1=m1, m2=m2)
  #options(warn = 0)
  
  LR.star = c(LR.obs, pmax(L, 0))
  
  # FIXME: make this optional or relocate it
  #print(summary(LR.star))
  return(LR.star)
  
}

simulateone = function(m) {

  iboot <- sample.int(m$n, replace=TRUE)
  
  # response
  #y = m$obj$env$data$yobs
  # for some reason yobs has NaNs and values differ from supplied data
  # FIXME is this a bug ?
  
  # in the meantime use frame
  respCol <- attr(terms(m$frame), "response")
  y <- m$frame[,respCol]
  n <- m$n
  
  # get mus: without ranefs from report, with ranefs, simulate  
  # simulate ranefs (parametric bootstrap)
  mus <- m$obj$env$simulate()$mu_predict
  
  # resampling rows function
  resamp = function(pits) {
    as.vector(matrix(pits, nrow=n)[iboot,])
  }
  
  # get dispersion parameters if any
  # FIXME, should check if any dispersion info rather than look at names
  #if (any(names(m$fit$par)=="betad")) {
  phi <- with(m$obj$env, exp(data$Xd %*% parList()$betad + ifelse(is.null(data$doffset), 0, data$doffset)))
  #}
  y.star <- switch(m$modelInfo$family$family, 
                  "gaussian" = {qnorm(resamp(pnorm(y, mus, phi)), mus, phi)},
                  "poisson" = {w = runif(n); qpois(resamp(pmin(w*ppois(y-1, mus) + (1-w)*ppois(y, mus), 1-1e-8)), mus)},
                  "nbinom1" = {0},
                  "nbinom2" = {w = runif(n); qnbinom(resamp(pmin(w * pnbinom(y-1, mu = mus, size = phi) + (1-w) * pnbinom(y, mu = mus, size = phi), 1-1e-8)), mu = mus, size=phi)},
                  "binomial" = {rbinom(resamp(runif(length(y), min = pbinom(y-1, 1, mus), max = pbinom(y, 1, mus))), 1, mus)})
  
  if (!all(is.finite(y.star)))
    stop("Infinite values in PITtrap")
  as.numeric(y.star)       
  
}

# refit
refitone = function(m1, m2, y.star) {
  
  # change data in TMBstruct 
  m1$TMBstruct$data.tmb$yobs = y.star
  m2$TMBstruct$data.tmb$yobs = y.star
  
  #m1$obj$env$data$yobs = y.star
  #m2$obj$env$data$yobs = y.star
  
  #m1$obj$par = with(m1$obj$env, last.par.best[-random])
  #m2$obj$par = with(m2$obj$env, last.par.best[-random])
  
  # weaken convergence criteria
  m1$TMBstruct$control$optCtrl[["rel.tol"]] = 1e-3
  m2$TMBstruct$control$optCtrl[["rel.tol"]] = 1e-3
  
  # warm start from bestfit parameters
  m1$TMBstruct$parameters = m1$obj$env$parList()
  m2$TMBstruct$parameters = m2$obj$env$parList()
  
  # possibly change control parameters to speed thinks up
  # refit
  m1.star = glmmTMB:::fitTMB(m1$TMBstruct)
  m2.star = glmmTMB:::fitTMB(m2$TMBstruct)
  #fit1 = do.call("optim", c(m1$obj, list(control=list(abstol=1e-3, reltol=1e-3))))
  #fit2 = do.call("optim", c(m2$obj, list(control=list(abstol=1e-3, reltol=1e-3))))  
  # return log of likelihood ratio
  #return(fit1$value - fit2$value)
  return(m1.star$fit$objective - m2.star$fit$objective)
}

# residuals, with PIT-residuals added

##' Compute residuals for a nuglmm object
##'
##' @param object a \dQuote{nuglmm} object
##' @param type (character) residual type
##' @param \dots ignored, for method compatibility
##' @importFrom stats fitted predict model.response residuals
##' @export
residuals.nuglmm <- function(object, type=c("response", "pearson", "PIT"), ...) {
  type <- match.arg(type)
  if (type == "PIT") {

    # get response
    y <- model.response(object$frame)

    # get mus
    mus <- predict(object,type="response")
    
    n <- object$n
    
    # dispersion
    phi <- with(object$obj$env, exp(data$Xd %*% parList()$betad + ifelse(is.null(data$doffset), 0, data$doffset)))

    # try sum on log scale to improve numerical stability ?
    # see https://andrewgelman.com/2016/06/11/log-sum-of-exponentials/
    log_sum_exp <- function(u, v) {max(u, v) + log(exp(u - max(u, v)) + exp(v - max(u, v)))}
    squeeze <- function(u) {(1.0 - .Machine$double.eps) * (u - .5) + .5}
    
    switch(object$modelInfo$family$family, 
                    "gaussian" = {pnorm(y, mus, phi)},
                    "poisson" = {w = runif(n); pmin(w*ppois(y-1, mus) + (1-w)*ppois(y, mus), 1-1e-8)},
                    "nbinom2" = {w = runif(n); pmin(w * pnbinom(y-1, mu = mus, size = phi) + (1-w) * pnbinom(y, mu = mus, size = phi), 1-1e-8)},
                    "binomial" = {w = runif(n); w * pbinom(y-1, 1, mus) + (1-w) * pbinom(y, 1, mus)})
    
  } else {
    matrix(glmmTMB:::residuals.glmmTMB(object, type, ...), nrow = object$n)
  }
}