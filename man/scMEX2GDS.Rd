\name{scMEX2GDS}
\alias{scMEX2GDS}
\title{Convert MEX files to GDS}
\description{
    Creates a single-cell GDS file from Cell Ranger MEX files.
}
\usage{
scMEX2GDS(feature_fn, barcode_fn, mtx_fn, outfn,
    feature_colnm=c("id", "gene", "feature_type"),
    type=c("float32", "float64", "int32"), compress="LZMA_RA", clean=TRUE,
    verbose=TRUE)
}
\arguments{
    \item{feature_fn}{the input file name for features}
    \item{barcode_fn}{the input file name for barcodes}
    \item{mtx_fn}{the input count matrix in MEX format}
    \item{outfn}{the output file name in GDS format}
    \item{feature_colnm}{the column names used in \code{feature_fn}}
    \item{type}{numeric data type in the output file}
    \item{compress}{the compression method, see \code{\link{add.gdsn}}}
    \item{clean}{TRUE}
    \item{verbose}{if \code{TRUE}, show information}
}
\value{
    Return the path of the output file.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{scConvGDS}}, \code{\link{scHDF2GDS}}
}

\keyword{CellRanger}
\keyword{GDS}
\keyword{SingleCell}
