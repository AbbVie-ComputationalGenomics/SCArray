\name{scExperiment}
\alias{scExperiment}
\title{Get a SummarizedExperiment}
\description{
    Gets an instance of SingleCellExperiment or SummarizedExperiment.
}
\usage{
scExperiment(gdsfile, sce=TRUE, use.names=TRUE, load.row=TRUE, load.col=TRUE)
}
\arguments{
    \item{gdsfile}{character for a file name, or a single-cell GDS object with
        class \code{SCArrayFileClass}}
    \item{sce}{if \code{TRUE}, return an instance of
        \code{SingleCellExperiment}, otherwise an instance of
        \code{SummarizedExperiment}}
    \item{use.names}{if \code{TRUE}, load \code{dimnames} from 'feature.id' and
        'sample.id'}
    \item{load.row}{\code{TRUE} for loading \code{rowData} from the gds node
        "feature.data" in \code{gdsfile}}
    \item{load.col}{\code{TRUE} for loading \code{colData} from the gds node
        "sample.data" in \code{gdsfile}}
}
\value{
    Return an instance of \code{\link{SingleCellExperiment}} or
\code{\link{SummarizedExperiment}}.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{scOpen}}, \code{\link{scClose}}
}

\examples{
# a GDS file for SingleCellExperiment
fn <- system.file("extdata", "example.gds", package="SCArray")

sce <- scExperiment(fn)
sce

remove(sce)
}

\keyword{GDS}
\keyword{SingleCell}
