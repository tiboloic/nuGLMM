\name{anova.nuglmm}
\alias{anova.nuglmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
anova.nuglmm(object, ..., nboot = 0, model.names = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
  \item{nboot}{
%%     ~~Describe \code{nboot} here~~
}
  \item{model.names}{
%%     ~~Describe \code{model.names} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (object, ..., nboot = 0, model.names = NULL) 
{
    mCall <- match.call(expand.dots = TRUE)
    dots <- list(...)
    .sapply <- function(L, FUN, ...) unlist(lapply(L, FUN, ...))
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
        if (nboot > 0) {
            null.dists <- matrix(0, nboot + 1, length(mods) - 
                1)
            for (imod in 1:(length(mods) - 1)) null.dists[, imod] <- pittrap(mods[[imod]], 
                mods[[imod + 1]], nboot, ...)
            val <- data.frame(Df = Df, AIC = .sapply(llks, AIC), 
                BIC = .sapply(llks, BIC), logLik = llk, deviance = -2 * 
                  llk, LR = c(NA, null.dists[1, ]), `Pr(>LR)` = c(NA, 
                  apply(null.dists, 2, function(null.dist) mean(null.dist[1] <= 
                    null.dist))), row.names = names(mods), check.names = FALSE)
            header <- c(paste0("Null distribution of likelihood ratio estimated by ", 
                nboot, " PITtrap simulations"), header)
        }
        else {
            val <- data.frame(Df = Df, AIC = .sapply(llks, AIC), 
                BIC = .sapply(llks, BIC), logLik = llk, deviance = -2 * 
                  llk, Chisq = chisq, `Chi Df` = dfChisq, `Pr(>Chisq)` = pchisq(chisq, 
                  dfChisq, lower.tail = FALSE), row.names = names(mods), 
                check.names = FALSE)
        }
        class(val) <- c("anova", class(val))
        forms <- lapply(lapply(calls, `[[`, "formula"), deparse)
        ziforms <- lapply(lapply(calls, `[[`, "ziformula"), deparse)
        dispforms <- lapply(lapply(calls, `[[`, "dispformula"), 
            deparse)
        structure(val, heading = c(header, "Models:", paste(paste(paste(rep(names(mods), 
            times = lengths(forms)), unlist(forms), sep = ": "), 
            unlist(ziforms), sep = ", zi="), unlist(dispforms), 
            sep = ", disp=")))
    }
    else stop("no single-model anova() method for nuglmm")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
