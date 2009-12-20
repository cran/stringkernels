\name{gapweightkernel}
\Rdversion{1.1}
\alias{gapweightkernel}
\title{
  Gap-weighted string kernel
}
\description{
  Creates a kernel object for the gap-weighted string kernel. 
  This kernel uses words as tokens by default. Use in conjunction with \pkg{kernlab}.
}
\usage{
gapweightkernel(length = 2, lambda = 0.75, normalized = TRUE, 
    tokenizer = openNLP::tokenize, use_characters = FALSE)
}
\arguments{
  \item{length}{
    Match length (excluding gaps)
  }
  \item{lambda}{
    Gap length penalty factor
  }
  \item{normalized}{
    Normalize kernel values (default: \code{TRUE})
  }
  \item{tokenizer}{
    String tokenizer function. By default, this uses \pkg{openNLP}'s \code{tokenize} to split 
    the text into words, but users may specify their own function.
    Ignored if \code{use_characters} is \code{TRUE}.
  }
  \item{use_characters}{
    Split texts by character, rather than by word.
  }
}
\details{
This kernel generation function returns a kernel that computes the number of 
gapped (non-contiguous) matches of \code{length} matching tokens between two strings. Gaps are penalized by a factor 
\code{lambda}, i.e., each match is assigned a weight of  
\eqn{\lambda^{L - l}}{lambda ^ (L - l)}, with \eqn{L} as the total match length and \eqn{l} as \code{length}.

By default, this kernel uses words (and punctuation marks) rather than characters as atomic tokens. 
This is usually yields better results than gapped character matching.

This implementation is based on the gapped substring kernel by Rousu/Shawe-Taylor. 
Note that this algorithm is optimized for large alphabets, usually consisting of words. 

}
\value{
An S4 object of class \code{stringKernelEx}.
}
\references{
Juho Rousu and John Shawe-Taylor. Efficient computation of gapped substring
kernels on large alphabets. Journal of Machine Learning Research,
6:1323-1344, 2005.
}
\author{
Martin Kober \cr \email{martin.kober@gmail.com}
}


\seealso{
\code{\link{multigapweightkernel}}
}
\examples{

s = "The cat was chased by the fat dog"
t = "The fat cat bit the dog"
gwk = gapweightkernel()
gwk(s,t)

gwk2 = gapweightkernel(length=4, normalized=FALSE)
gwk2(s,t)

gwk3 = gapweightkernel(lambda=1, normalized=FALSE)
gwk3(s,t)


}
\keyword{misc}
