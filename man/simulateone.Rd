\name{simulateone}
\alias{simulateone}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simulateone(m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{
%%     ~~Describe \code{m} here~~
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
function (m) 
{
    iboot <- sample.int(m$n, replace = TRUE)
    respCol <- attr(terms(m$frame), "response")
    y <- m$frame[, respCol]
    n <- m$n
    mus <- m$obj$env$simulate()$mu_predict
    resamp = function(pits) {
        as.vector(matrix(pits, nrow = n)[iboot, ])
    }
    phi <- with(m$obj$env, exp(data$Xd \%*\% parList()$betad + 
        ifelse(is.null(data$doffset), 0, data$doffset)))
    y.star <- switch(m$modelInfo$family$family, gaussian = {
        qnorm(resamp(pnorm(y, mus, phi)), mus, phi)
    }, poisson = {
        w = runif(n)
        qpois(resamp(pmin(w * ppois(y - 1, mus) + (1 - w) * ppois(y, 
            mus), 1 - 1e-08)), mus)
    }, nbinom1 = {
        0
    }, nbinom2 = {
        w = runif(n)
        qnbinom(resamp(pmin(w * pnbinom(y - 1, mu = mus, size = phi) + 
            (1 - w) * pnbinom(y, mu = mus, size = phi), 1 - 1e-08)), 
            mu = mus, size = phi)
    }, binomial = {
        rbinom(resamp(runif(length(y), min = pbinom(y - 1, 1, 
            mus), max = pbinom(y, 1, mus))), 1, mus)
    })
    if (!all(is.finite(y.star))) 
        stop("Infinite values in PITtrap")
    as.numeric(y.star)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
