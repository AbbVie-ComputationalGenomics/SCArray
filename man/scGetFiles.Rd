\name{scGetFiles}
\alias{scGetFiles}
\alias{scGetFiles,SC_GDSArray-method}
\alias{scGetFiles,SummarizedExperiment-method}

\title{File names for on-disk backend}
\description{
    Get a list of file names for DelayedArray with an on-disk backend.
}

\usage{
scGetFiles(object, ...)
\S4method{scGetFiles}{SC_GDSArray}(object, ...)
\S4method{scGetFiles}{SummarizedExperiment}(object, ...)
}

\arguments{
    \item{object}{input R object (e.g., a GDS-specific DelayedArray)}
    \item{...}{additional arguments passed to specific methods}
}
\value{
    Return a character vector storing file names.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{path}}
}

\keyword{methods}
\keyword{GDS}
