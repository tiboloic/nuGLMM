\name{residuals.nuglmm}
\alias{residuals.nuglmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
residuals.nuglmm(object, type = c("response", "pearson", "PIT"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
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
function (object, type = c("response", "pearson", "PIT"), ...) 
{
    type <- match.arg(type)
    if (type == "PIT") {
        y <- model.response(object$frame)
        mus <- predict(object, type = "response")
        n <- object$n
        phi <- with(object$obj$env, exp(data$Xd \%*\% parList()$betad + 
            ifelse(is.null(data$doffset), 0, data$doffset)))
        log_sum_exp <- function(u, v) {
            max(u, v) + log(exp(u - max(u, v)) + exp(v - max(u, 
                v)))
        }
        squeeze <- function(u) {
            (1 - .Machine$double.eps) * (u - 0.5) + 0.5
        }
        switch(object$modelInfo$family$family, gaussian = {
            pnorm(y, mus, phi)
        }, poisson = {
            w = runif(n)
            pmin(w * ppois(y - 1, mus) + (1 - w) * ppois(y, mus), 
                1 - 1e-08)
        }, nbinom2 = {
            w = runif(n)
            pmin(w * pnbinom(y - 1, mu = mus, size = phi) + (1 - 
                w) * pnbinom(y, mu = mus, size = phi), 1 - 1e-08)
        }, binomial = {
            w = runif(n)
            w * pbinom(y - 1, 1, mus) + (1 - w) * pbinom(y, 1, 
                mus)
        })
    }
    else {
        matrix(glmmTMB:::residuals.glmmTMB(object, type, ...), 
            nrow = object$n)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
