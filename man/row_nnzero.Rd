\name{row_nnzero}
\alias{row_nnzero}
\alias{row_nnzero,matrix-method}
\alias{row_nnzero,Matrix-method}
\alias{row_nnzero,DelayedMatrix-method}
\alias{row_nnzero,SC_GDSMatrix-method}
\alias{col_nnzero}
\alias{col_nnzero,matrix-method}
\alias{col_nnzero,Matrix-method}
\alias{col_nnzero,DelayedMatrix-method}
\alias{col_nnzero,SC_GDSMatrix-method}

\title{Numbers of Non-zeros}
\description{
    Calculates the numbers of non-zeros for each row or column of a matrix-like
object.
}

\usage{
row_nnzero(x, na.counted=NA, ...)
col_nnzero(x, na.counted=NA, ...)

\S4method{row_nnzero}{matrix}(x, na.counted=NA, ...)
\S4method{row_nnzero}{Matrix}(x, na.counted=NA, ...)
\S4method{row_nnzero}{DelayedMatrix}(x, na.counted=NA, ...)
\S4method{row_nnzero}{SC_GDSMatrix}(x, na.counted=NA, ...)

\S4method{col_nnzero}{matrix}(x, na.counted=NA, ...)
\S4method{col_nnzero}{Matrix}(x, na.counted=NA, ...)
\S4method{col_nnzero}{DelayedMatrix}(x, na.counted=NA, ...)
\S4method{col_nnzero}{SC_GDSMatrix}(x, na.counted=NA, ...)
}

\arguments{
    \item{x}{a matrix-like object}
    \item{na.counted}{a logical: TRUE for counting NA/NaN as non-zero, FALSE
        for counting NA/NaN as zero, NA (default) for return NA when
        encountering NA/NaN}
    \item{...}{additional arguments passed to specific methods}
}
\value{
    Return an integer vector object for the numbers of non-zeros.
}

\author{Xiuwen Zheng}

\examples{
# a GDS file for SingleCellExperiment
fn <- system.file("extdata", "example.gds", package="SCArray")

cnt <- scArray(fn, "counts")
cnt

row_nnzero(cnt, na.counted=TRUE)
col_nnzero(cnt, na.counted=TRUE)


rm(cnt)
}

\keyword{methods}
\keyword{GDS}
