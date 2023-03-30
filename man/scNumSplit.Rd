\name{scNumSplit}
\alias{scNumSplit}
\title{Split a number}
\description{
    Splits a number into multiple groups with equal size.
}
\usage{
scNumSplit(num, BPPARAM=getAutoBPPARAM())
}
\arguments{
    \item{num}{a length-one number (the total count) for splitting
        (must be >= 0)}
    \item{BPPARAM}{NULL, a number for the number of groups, or a
        \code{BiocParallelParam} object; if not specified, call
        \code{getAutoBPPARAM()}}
}
\value{
    Return a list of length-two numeric vectors for the start and end
positions. \code{BPPARAM=NULL} is as the same as \code{BPPARAM=1}, if it is a
\code{BiocParallelParam} object, call \code{bpnworkers()} to get the number of
groups.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{getAutoBPPARAM}}, \code{\link{BiocParallelParam}},
    \code{\link{bpnworkers}}
}

\examples{
scNumSplit(100, NULL)
scNumSplit(100, 0)
scNumSplit(100, 1)
scNumSplit(100, 3)
scNumSplit(100)

scNumSplit(0)  # zero-length
}

\keyword{GDS}
