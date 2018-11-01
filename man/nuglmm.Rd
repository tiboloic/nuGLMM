\name{nuglmm}
\alias{nuglmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
nuglmm(formula, data = NULL, offset = NULL, weights = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{offset}{
%%     ~~Describe \code{offset} here~~
}
  \item{weights}{
%%     ~~Describe \code{weights} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (formula, data = NULL, offset = NULL, weights = NULL, 
    ...) 
{
    cl <- match.call()
    environment(formula) <- parent.frame()
    lhs <- eval(formula[[2]], data, enclos = environment(formula))
    if (!is.null(lhs)) {
        if (!is.matrix(lhs)) {
            stop("nuglmm needs a matrix-valued response")
        }
        else {
            if (ncol(lhs) < 5) 
                stop("A minimumn of 5 columns (species) is required")
        }
    }
    else stop("Cannot evaluate response")
    p <- ncol(lhs)
    n <- nrow(lhs)
    f <- formula(Y[, 1] ~ 1)
    f[[2]][[2]] <- formula[[2]]
    f[[3]] <- formula[[3]]
    m <- tryCatch(glmmTMB(f, data, doFit = FALSE, ...), error = function(e) e)
    if ("error" \%in\% class(m)) 
        stop(m$message)
    .valid.family <- c("gaussian", "poisson", "nbinom2", "binomial")
    if (!m$family$family \%in\% .valid.family) 
        stop(paste0("This implementation of nuglmm accepts only the following families: ", 
            .valid.family))
    fr <- m$fr
    respCol <- attr(terms(fr), "response")
    names(respCol) <- deparse(formula[[2]])
    y <- lhs
    fr. = as.data.frame(lapply(fr[-respCol], rep, p))
    fr.[names(respCol)] = as.vector(y)
    fr.[".spp"] = gl(p, nrow(fr))
    nuformula = formula
    makeOp <- function(x, y, op = NULL) {
        if (is.null(op) || missing(y)) {
            if (is.null(op)) {
                substitute(OP(X), list(X = x, OP = y))
            }
            else {
                substitute(OP(X), list(X = x, OP = op))
            }
        }
        else substitute(OP(X, Y), list(X = x, OP = op, Y = y))
    }
    newTerms <- list(formula)
    reTerms = lme4:::findbars(formula)
    for (ireTerm in seq_along(reTerms)) {
        reTerm = reTerms[[ireTerm]]
        newTerms[[ireTerm + 1]] = makeOp(makeOp(reTerm[[2]], 
            makeOp(reTerm[[3]], quote(.spp), quote(`:`)), reTerm[[1]]), 
            quote(`(`))
    }
    newTerms[[length(newTerms) + 1]] <- makeOp(makeOp(lme4:::nobars(formula)[[3]], 
        quote(.spp), quote(`|`)), op = quote(diag))
    nuformula <- Reduce(function(f, term) {
        glmmTMB:::RHSForm(f) <- makeOp(glmmTMB:::RHSForm(f), 
            term, quote(`+`))
        return(f)
    }, newTerms)
    if (!missing(offset)) {
        if (is.matrix(offset)) {
            if (all(dim(offset) == c(n, p))) 
                offset <- as.vector(offset)
            else stop("Incorrect dimension of offset matrix")
        }
    }
    if (!missing(weights)) {
        if (is.matrix(weights)) {
            if (all(dim(weights) == c(n, p))) 
                weights <- as.vector(weights)
            else stop("Incorrect dimension of weight matrix")
        }
    }
    m.struct <- glmmTMB(nuformula, data = fr., offset = offset, 
        weights = weights, doFit = FALSE, ...)
    m.struct$data.tmb$whichPredict <- 1:nrow(fr.)
    m <- glmmTMB:::fitTMB(m.struct)
    m$n <- nrow(fr)
    m$p <- p
    m$TMBstruct <- m.struct
    class(m) <- c("nuglmm", class(m))
    m
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
