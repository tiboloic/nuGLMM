\name{pittrap}
\alias{pittrap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pittrap(m1, m2, nboot, ..., parallel = c("no", "multicore", "snow"), ncpus = getOption("nuglmm.ncpus"), cl = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m1}{
%%     ~~Describe \code{m1} here~~
}
  \item{m2}{
%%     ~~Describe \code{m2} here~~
}
  \item{nboot}{
%%     ~~Describe \code{nboot} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
  \item{parallel}{
%%     ~~Describe \code{parallel} here~~
}
  \item{ncpus}{
%%     ~~Describe \code{ncpus} here~~
}
  \item{cl}{
%%     ~~Describe \code{cl} here~~
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
function (m1, m2, nboot, ..., parallel = c("no", "multicore", 
    "snow"), ncpus = getOption("nuglmm.ncpus"), cl = NULL) 
{
    if (missing(ncpus)) 
        ncpus <- parallel::detectCores()
    if (missing(parallel)) {
        parallel <- getOption("nuglmm.parallel")
        if (is.null(parallel)) {
            parallel <- ifelse(parallel::detectCores() > 1, "multicore", 
                "no")
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
    }
    LR.obs <- m1$fit$objective - m2$fit$objective
    m1$TMBstruct$parameters <- m1$obj$env$parList()
    m2$TMBstruct$parameters <- m2$obj$env$parList()
    if (have_snow) 
        FUN = function(i, m1, m2) {
            require(glmmTMB)
            source("Methods.R")
            y.star = simulateone(m1)
            return(refitone(m1, m2, y.star))
        }
    else FUN = function(i) {
        y.star = simulateone(m1)
        return(refitone(m1, m2, y.star))
    }
    if (ncpus > 1) {
        if (have_mc) {
            L <- parallel::mclapply(1:nboot, FUN, mc.cores = ncpus)
        }
        else if (have_snow) {
            if (is.null(cl)) {
                new_cl <- TRUE
                cl <- parallel::makePSOCKcluster(rep("localhost", 
                  ncpus))
            }
            L <- parallel::clusterApply(cl, 1:nboot, FUN, m1 = m1, 
                m2 = m2)
            if (new_cl) {
                parallel::stopCluster(cl)
            }
        }
    }
    else {
        L <- lapply(1:nboot, FUN)
    }
    L <- unlist(L)
    LR.star = c(LR.obs, pmax(L, 0))
    return(LR.star)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
