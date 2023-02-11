\name{scSetBounds}
\alias{scSetBounds}
\alias{scSetMax}
\alias{scSetMin}
\title{Set the bounds}
\description{
    Set the maximum and/or minimum on a GDS-specific DelayedArray.
}
\usage{
scSetMax(x, vmax)
scSetMin(x, vmin)
scSetBounds(x, vmin=NaN, vmax=NaN)
}
\arguments{
    \item{x}{a \code{SC_GDSArray} object}
    \item{vmax}{maximum, length-one}
    \item{vmin}{minimum, length-one}
}
\value{
    Return an object with the class \code{SC_GDSMatrix} or \code{SC_GDSArray}.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{scReplaceNA}}
}

\examples{
suppressPackageStartupMessages(library(DelayedArray))

m <- matrix(1:12, nrow=3)
(mat <- DelayedArray(m))

new_m <-  scObj(mat)  # wrap a in-memory DelayedMatrix
class(new_m)  # SC_GDSMatrix

scSetMax(new_m, 5)
scSetMin(new_m, 5)
scSetBounds(new_m, 4, 9)
}

\keyword{GDS}
