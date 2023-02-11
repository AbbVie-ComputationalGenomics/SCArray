\name{scHDF2GDS}
\alias{scHDF2GDS}
\title{Convert HDF5 files to GDS}
\description{
    Creates a single-cell GDS file from Cell Ranger HDF5 files.
}
\usage{
scHDF2GDS(h5_fn, outfn, group=c("matrix", "mm10"), feature_path=character(),
    type=c("float32", "float64", "int32"), compress="LZMA_RA", clean=TRUE,
    verbose=TRUE)
}
\arguments{
    \item{h5_fn}{the input HDF5 file name}
    \item{outfn}{the output file name in GDS format}
    \item{group}{the name of the group in the HDF5 file where the sparse
        matrix is stored; if there are more than one group names, the first
        existing group in the HDF5 file is used; \code{"mm10"} is usually
        used for 10x Genomics datasets}
    \item{feature_path}{a character vector for feature variables, otherwise
        detecting automatically}
    \item{type}{numeric data type in the output file}
    \item{compress}{the compression method, see \code{\link{add.gdsn}}}
    \item{clean}{TRUE}
    \item{verbose}{if \code{TRUE}, show information}
}
\details{
    The packages \pkg{rhdf5} and \pkg{HDF5Array} should be installed.
}
\value{
    Return the path of the output file.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{scConvGDS}}, \code{\link{scMEX2GDS}}
}

\keyword{CellRanger}
\keyword{GDS}
\keyword{SingleCell}
