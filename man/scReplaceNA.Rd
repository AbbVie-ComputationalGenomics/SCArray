\name{scReplaceNA}
\alias{scReplaceNA}
\title{Replacement}
\description{
    Replace NA/NaN in a GDS-specific DelayedArray by a specified value.
}
\usage{
scReplaceNA(x, v=0L)
}
\arguments{
    \item{x}{a \code{SC_GDSArray} object}
    \item{v}{a length-one double or integer value}
}
\value{
    Return an object with the class \code{SC_GDSMatrix} or \code{SC_GDSArray}.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{scSetMin}}, \code{\link{scSetMax}}, \code{\link{scSetBounds}}
}

\examples{
suppressPackageStartupMessages(library(DelayedArray))

m <- matrix(1:12, nrow=3)
m[2, c(1,3)] <- NA
(mat <- DelayedArray(m))

new_m <-  scObj(mat)  # wrap a in-memory DelayedMatrix
class(new_m)  # SC_GDSMatrix

scReplaceNA(new_m,  999)
}

\keyword{GDS}
