\name{SCArray-package}
\alias{SCArray-package}
\alias{SCArray}
\docType{package}
\title{Large-scale single-cell RNA-seq data manipulation with GDS files}
\description{
    The package combines dense/sparse matrices stored in GDS files and the
Bioconductor infrastructure framework to provide out-of-memory data storage
and manipulation using the R programming language.
}
\details{
    \tabular{ll}{
        Package: \tab SCArray \cr
        Type: \tab Package \cr
        License: \tab GPL version 3 \cr
    }
}
\author{
    Xiuwen Zheng \email{xiuwen.zheng@abbvie.com}
}

\examples{
# a GDS file for SingleCellExperiment
fn <- system.file("extdata", "example.gds", package="SCArray")

sce <- scExperiment(fn)
sce

rm(sce)
}

\keyword{GDS}
\keyword{SingleCell}
