\name{multigapweightkernel}
\Rdversion{1.1}
\alias{multigapweightkernel}
\alias{multigapweight-class}
\alias{getkernel}
\alias{getkernel-methods}
\alias{getkernel,multigapweight-method}

\title{
Multiple gap-weight kernels
}
\description{
Compute gap-weight kernels of multiple length at once and pack them in a precomputed kernel.
}
\usage{
multigapweightkernel(items, maxlength, kernelarray = NULL, lambda = 0.75, 
    normalized = TRUE, tokenizer = openNLP::tokenize, minlength = 1)

\S4method{getkernel}{multigapweight}(mgw, length, use_dummy = FALSE)
}
\arguments{
  \item{items}{
    List of input texts
}
  \item{maxlength}{
    Maximum match length
}
  \item{kernelarray}{
    Optionally supply an array of kernel values.
}
  \item{lambda}{
    Gap length penalty factor
}
  \item{normalized}{
    Normalize kernel values
}
  \item{tokenizer}{
    String tokenizer function. By default, this uses \pkg{openNLP}'s \code{tokenize} to split 
    the text into words, but users may specify their own function.
}
  \item{minlength}{
    Minimum match length
}
  \item{mgw}{
    multigapweight object returned by \code{multigapweightkernel}
}
  \item{length}{
    The desired length parameter for the kernel 
}
  \item{use_dummy}{
    The flag \code{use_dummy=TRUE} can be used to 
    create a kernel with dummy values (see \code{\link{precomputedkernel}})
}
}
\details{
The dynamic programming algorithm used for the gap-weighted kernel works by
computing the matching statistics for an incrementally larger match length.
 
Therefore, computing the kernel value for match length \eqn{n} does not take 
significantly less computational time than computing all kernel values for 
\eqn{n'\leq n}{n' <= n}. 

This function computes kernel matrices for multiple lengths in one step.
The \code{getkernel} method retrieves the matrix of the desired length and 
creates a kernel object with the precomputed values.

}
\value{
A \code{multigapweight} object that contains the kernel value array (a kernel matrix with
an additional dimension for length) and the kernel parameters.
}

\author{
Martin Kober \cr \email{martin.kober@gmail.com}
}

\seealso{
\code{\link{precomputedkernel}}
}
\examples{

library(tm)
data(crude)

m = multigapweightkernel(crude, maxlength=3, minlength=2)

k2 = getkernel(m, 2)
k3 = getkernel(m, 3)

kernelMatrix(k2, crude[1:5])
kernelMatrix(k3, crude[1:5])


}
\keyword{ misc }
