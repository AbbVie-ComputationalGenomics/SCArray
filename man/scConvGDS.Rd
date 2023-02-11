\name{scConvGDS}
\alias{scConvGDS}
\title{Create a GDS file}
\description{
    Creates a single-cell GDS file from an R object.
}
\usage{
scConvGDS(obj, outfn, assay.name=NULL, save.sp=TRUE,
    type=c("float32", "float64", "int32"), compress="LZMA_RA", clean=TRUE, verbose=TRUE)
}
\arguments{
    \item{obj}{a dense/sparse matrix, DelayedMatrix, SummarizedExperiment or
        SingleCellExperiment}
    \item{outfn}{the output file name in GDS format}
    \item{assay.name}{a character vector for assay names or \code{NULL}; if \code{NULL},
        to include all available assays, otherwise only include the assays in
        \code{assay.name}}
    \item{save.sp}{if \code{TRUE}, save it to a sparse matrix in GDS; otherwise,
        store dense matrix}
    \item{type}{numeric data type in the output file}
    \item{compress}{the compression method, see \code{\link{add.gdsn}}; or
        \code{""} for no data compression}
    \item{clean}{TRUE}
    \item{verbose}{if \code{TRUE}, show information}
}
\value{
    Return the path of the output file.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{scOpen}}, \code{\link{scClose}},
    \code{\link{scMEX2GDS}}, \code{\link{scHDF2GDS}}
}

\examples{
# load a SingleCellExperiment object
fn <- system.file("extdata", "example.rds", package="SCArray")
sce <- readRDS(fn)
sce

scConvGDS(sce, "test.gds")

# remove the temporary file
unlink("test.gds")
}

\keyword{GDS}
\keyword{SingleCell}
