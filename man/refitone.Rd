\name{refitone}
\alias{refitone}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
refitone(m1, m2, y.star)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m1}{
%%     ~~Describe \code{m1} here~~
}
  \item{m2}{
%%     ~~Describe \code{m2} here~~
}
  \item{y.star}{
%%     ~~Describe \code{y.star} here~~
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
function (m1, m2, y.star) 
{
    m1$TMBstruct$data.tmb$yobs = y.star
    m2$TMBstruct$data.tmb$yobs = y.star
    m1$TMBstruct$control$optCtrl[["rel.tol"]] = 0.001
    m2$TMBstruct$control$optCtrl[["rel.tol"]] = 0.001
    m1$TMBstruct$parameters = m1$obj$env$parList()
    m2$TMBstruct$parameters = m2$obj$env$parList()
    m1.star = glmmTMB:::fitTMB(m1$TMBstruct)
    m2.star = glmmTMB:::fitTMB(m2$TMBstruct)
    return(m1.star$fit$objective - m2.star$fit$objective)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
