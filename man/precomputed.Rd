\name{precomputedkernel}
\Rdversion{1.1}
\alias{precomputedkernel}
\alias{precomputeddummy}
\alias{stringkernelPrecomputed-class}
\alias{kernelMatrix,stringkernelPrecomputed-method}

\title{
Kernel with precomputed values.
}
\description{
This function creates a kernel that stores precomputed kernel values and 
retrieves them as needed.
Use in conjunction with \pkg{kernlab}.
}
\usage{
precomputedkernel(kernel, items, kernelmatrix = NULL, use_kernel = FALSE, 
    use_dummy = TRUE)

precomputeddummy(items)
}
\arguments{
  \item{kernel}{
    May be either a string kernel object of class \code{stringkernelEx} or an arbitrary 
    character string containing the name of the kernel used. In the latter case, 
    a \code{kernelmatrix} has to be supplied and \code{use_kernel} must be 
    \code{FALSE}.
}
  \item{items}{
    Items to be compared.
}
  \item{kernelmatrix}{
    A kernel matrix corresponding to the \code{items} and the \code{kernel}. 
    If not provided, \code{kernel} is used to compute it on-the-fly.
}
  \item{use_kernel}{
    If \code{TRUE}, \code{kernel} is called whenever the object is called with 
    texts not found in \code{items}. Otherwise, only stored values are used and 
    input values not in \code{items} will raise an error.
}
  \item{use_dummy}{
    Use dummy texts instead of the actual items. When true, the kernel does not 
    match the actual items, instead it interprets the input texts as integer indices. 
    See Details on how to use this. Cannot be used in conjunction with \code{use_kernel}.
}
}
\details{
On most string kernel tasks, computing the kernel is the most time-consuming operation.
This kernel can store kernel values and transparently return stored 
values instead of computing them on-the-fly. 
This is potentially useful whenever the same kernel values are needed
for multiple classification runs, e.g., cross-fold validation or parameter tuning.

The kernel internally matches input items to the stored items and retrieves 
the stored kernel values. When corpora are S4 objects (e.g., a corpus from \pkg{tm}) and/or very large, 
matching input items to stored items can consume a lot of time.
The \code{use_dummy} flag can be used to avoid this. The flag instructs the kernel function 
to interpret input values directly as indices into the kernel matrix.
These ``dummy'' input values can be created with \code{precomputeddummy(items)}.
The dummy values are compatible with all training and prediction functions in \pkg{kernlab}
that expect strings as input (see Examples).



}
\value{
A S4 kernel object of class \code{stringkernelPrecomputed}.
}
\author{
Martin Kober \cr \email{martin.kober@gmail.com}
}

\seealso{
\code{\link{multigapweightkernel}}
}
\examples{

library(tm)
library(kernlab)
data(crude)

wdk = worddot(type="spectrum", length=2)
kernelMatrix(wdk, crude[1:3], crude[17:20])

pre = precomputedkernel(wdk, crude, use_dummy=TRUE)
dummy = precomputeddummy(crude)

kernelMatrix(pre, dummy[1:3], dummy[17:20])

class = factor(rep(c(1,-1),10))
model = ksvm(dummy[1:10], class[1:10], kernel=pre)

predict(model, dummy[11:20])

}
\keyword{ misc }
