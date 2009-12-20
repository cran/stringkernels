\name{worddot}
\Rdversion{1.1}
\alias{worddot}
\title{
Word-based string kernels.
}
\description{
This function is analogous to \pkg{kernlab}'s \code{stringdot}, using words instead of characters.
}
\usage{
worddot(type = c("spectrum", "boundrange", "constant", "exponential"), 
    length = 4, lambda = 1.1, normalized = TRUE, 
    tokenizer = openNLP::tokenize)
}
\arguments{
  \item{type}{
    Type of kernel to be used. Four types are supported: \cr    
    \code{spectrum} Matches of exactly length \eqn{n}.\cr
    \code{boundrange} Matches of all lengths up to \eqn{n}. \cr
    \code{exponential} Matches of all lengths with exponentially decaying weighting. 
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
  \item{tokenizer}{
    String tokenizer function. By default, this uses \pkg{openNLP}'s \code{tokenize} to split 
    the text into words, but users may specify their own function.
}
}
\details{
  This function is identical to the \code{\link[kernlab]{stringdot}} function in \pkg{kernlab},
  only that it uses words instead characters as tokens.
}
\value{
An S4 kernel object of class \code{stringkernelEx}.
}
\author{
Martin Kober \cr \email{martin.kober@gmail.com}
}

\seealso{
\code{\link[kernlab]{stringdot}}
}
\examples{
s = "The cat was chased by the fat dog"
t = "The fat cat bit the dog"
wdk = worddot(type="spectrum", length=2, normalized=FALSE)
wdk(s,t)
wdk(s,s)

}

\keyword{misc}
