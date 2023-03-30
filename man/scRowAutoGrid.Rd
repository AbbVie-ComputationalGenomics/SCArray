\name{scRowAutoGrid}
\alias{scRowAutoGrid}
\alias{scColAutoGrid}
\title{Automatic grids for matrix-like objects}
\description{
    Create automatic grids (RegularArrayGrid or ArbitraryArrayGrid for sparse
matrices) to use for block processing of matrix-like objects, where the blocks
are made of full rows or full columns.
}
\usage{
scRowAutoGrid(x, force=FALSE, nnzero=NULL)
scColAutoGrid(x, force=FALSE, nnzero=NULL)
}
\arguments{
    \item{x}{a matrix-like object (e.g., a SC_GDSMatrix object)}
    \item{force}{a logical, only applicable when \code{x} is a sparse in-memory
        matrix or a sparse \code{SC_GDSMatrix} object, see details}
    \item{nnzero}{a numeric vector for the numbers of non-zeros for rows or
        columns, \code{NULL} (default) for calling \code{row_nnzero()}
        or \code{col_nnzero()} when needed}
}
\details{
    The functions return regular RegularArrayGrid (calling \code{rowAutoGrid()}
or \code{colAutoGrid}), when \code{x} is neither a sparse in-memory matrix nor
a sparse \code{SC_GDSMatrix} object; otherwise, make use of the information of
the numbers of non-zeros to create ArbitraryArrayGrid for more efficient grids.
    When \code{force} is applicable and \code{force=TRUE}, the functions
return \code{ArbitraryArrayGrid} which needs the \code{nnzero} values. For
\code{force=FALSE}, \code{scRowAutoGrid()} returns \code{ArbitraryArrayGrid}
when \code{x} is not transposed, and \code{scColAutoGrid()} returns
\code{ArbitraryArrayGrid} when \code{x} is transposed.
    If \code{nnzero=NULL} and it is needed, the numbers of non-zeros for rows
or columns will be calculated internally. For a large matrix, it is more
efficient when \code{nnzero} is pre-defined.
    The internal block size can be controlled by \code{setAutoBlockSize()}.
    If the number of blocks in \code{ArbitraryArrayGrid} is more than
\code{RegularArrayGrid}, the functions return \code{RegularArrayGrid} instead
when \code{force} is not \code{TRUE}.

    Usually, \code{gd <- scRowAutoGrid()} or \code{gd <- scColAutoGrid()} is
used together with \code{blockApply(, grid=gd, as.sparse=attr(gd, "as.sparse"))}
or \code{blockReduce(, grid=gd, as.sparse=attr(gd, "as.sparse"))} to take
advantage of sparse matrices.
}
\value{
    Return an object of \code{RegularArrayGrid} or \code{ArbitraryArrayGrid}.
\code{attr(, "as.sparse")} is a suggested logical value for \code{as.sparse}
in \code{blockApply()} or \code{blockReduce()}.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{rowAutoGrid}}, \code{\link{colAutoGrid}},
    \code{\link{setAutoBlockSize}},
    \code{\link{blockApply}}, \code{\link{blockReduce}}
}

\examples{
# a GDS file for SingleCellExperiment
fn <- system.file("extdata", "example.gds", package="SCArray")

cnt <- scArray(fn, "counts")
cnt

setAutoBlockSize(1048576)  # use 1MB

scRowAutoGrid(cnt)  # it returns RegularArrayGrid since cnt is not very sparse
rowAutoGrid(cnt)
scRowAutoGrid(cnt, force=TRUE)  # ArbitraryArrayGrid


library(Matrix)
cnt2 <- Diagonal(1e5)  # a very sparse matrix

scRowAutoGrid(cnt2)    # 5 blocks
length(rowAutoGrid(cnt2))    # 100000

scColAutoGrid(cnt2)          # 5 blocks
length(colAutoGrid(cnt2))    # 100000 blocks


setAutoBlockSize()  # reset

rm(cnt)
}

\keyword{GDS}
