\name{scMemory}
\alias{scMemory}
\alias{scMemory,DelayedArray-method}
\alias{scMemory,SummarizedExperiment-method}

\title{Load Data to Memory}
\description{
    Loads the internal data to memory for any on-disk object.
}

\usage{
scMemory(x, ...)
\S4method{scMemory}{DelayedArray}(x, ...)
\S4method{scMemory}{SummarizedExperiment}(x, ...)
}

\arguments{
    \item{x}{input R object (e.g., a DelayedArray)}
    \item{...}{additional arguments passed to specific methods}
}
\value{
    Return an object (it maybe a different type compared with \code{x}).
}

\author{Xiuwen Zheng}

\examples{
suppressPackageStartupMessages(library(DelayedArray))

m <- matrix(1:12, nrow=3)
(mat <- DelayedArray(m))

str(scMemory(mat))
}

\keyword{methods}
\keyword{GDS}
