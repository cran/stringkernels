\name{stringdotEx}
\Rdversion{1.1}
\alias{stringdotEx}
\title{
Character-based string kernels
}
\description{
Character-based string kernels, analogous to \pkg{kernlab}'s \code{stringdot}.
}
\usage{
stringdotEx(type = c("spectrum", "boundrange", "constant", "exponential"), 
    length = 4, lambda = 1.1, normalized = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{type}{
    Type of kernel to be used. Four types are supported: \cr
    \code{spectrum} Matches of exactly length \eqn{n}.\cr
    \code{boundrange} Matches of all lengths up to \eqn{n} \cr
    \code{exponential} Matches of all lengths with exponentially decaying weighting 
    \eqn{\lambda^{-n}}{lambda ^ (-n)}. \cr
    \code{constant} Matches of all lengths with equal weighting.

}
  \item{length}{
    Length of the substrings (only for spectrum and boundrange kernels)
}
  \item{lambda}{
    Weighting factor, must be \eqn{> 1} (only for exponential kernel)
}
  \item{normalized}{
    Normalize word kernel values (default: \code{TRUE})
}
}
\details{
  This function is identical to the \code{\link[kernlab]{stringdot}} function in \pkg{kernlab}, 
  but compatible with functions like \code{precomputedkernel} requiring a \code{stringkernelEx} object.
}
\value{
An S4 kernel object of class \code{stringkernelEx}.
}
\author{
Martin Kober \cr \email{martin.kober@gmail.com}
}
\note{
  All newline in the input strings are internally converted to carriage returns 
  (\code{\\r}). 
}

\seealso{
\code{\link[kernlab]{stringdot}}
}


\examples{
s = "The cat was chased by the fat dog"
t = "The fat cat bit the dog"
stk = stringdotEx(type="spectrum", length=3, normalized=FALSE)
stk(s,t)
stk(s,s)


}

\keyword{ misc }
