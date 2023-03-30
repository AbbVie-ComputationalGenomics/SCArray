\name{SCArray-utils}
\alias{SCArray-utils}

\alias{[}
\alias{[,SC_GDSArray-method}
\alias{[,SC_GDSArray,ANY,ANY,ANY-method}
\alias{[[}
\alias{[[,SC_GDSArray-method}
\alias{[[,SC_GDSArray,ANY,ANY-method}
\alias{Ops}
\alias{Ops,SC_GDSArray-method}
\alias{Math}
\alias{Math,SC_GDSArray-method}
\alias{+}
\alias{+,SC_GDSArray,missing-method}
\alias{-}
\alias{-,SC_GDSArray,missing-method}
\alias{names<-}
\alias{names<-,SC_GDSArray-method}
\alias{dimnames<-}
\alias{dimnames<-,SC_GDSArray,ANY-method}

\alias{scale}
\alias{scale,SC_GDSMatrix-method}
\alias{pmin2}
\alias{pmin2,SC_GDSArray,SC_GDSArray-method}
\alias{pmin2,SC_GDSArray,vector-method}
\alias{pmin2,vector,SC_GDSArray-method}
\alias{pmax2}
\alias{pmax2,SC_GDSArray,SC_GDSArray-method}
\alias{pmax2,SC_GDSArray,vector-method}
\alias{pmax2,vector,SC_GDSArray-method}

\title{SC_GDSArray subsetting, Ops, Math}
\description{
    Subsetting, Arith, Compare, Logic and Math operations on the SC_GDSArray
object.
}

\usage{
# x[i, j, ... , drop = TRUE]
\S4method{[}{SC_GDSArray}(i, j, ... , drop=TRUE)
# x[[i, j, ...]]
\S4method{[[}{SC_GDSArray}(i, j, ...)

\S4method{Ops}{SC_GDSArray}(e1, e2)
\S4method{Math}{SC_GDSArray}(x)

# names(x) <- value
# dimnames(x) <- value

# Centers and/or scales the columns of a matrix
\S4method{scale}{SC_GDSMatrix}(x, center=TRUE, scale=TRUE)

\S4method{pmin2}{SC_GDSArray,SC_GDSArray}(e1, e2)
\S4method{pmin2}{SC_GDSArray,vector}(e1, e2)
\S4method{pmin2}{vector,SC_GDSArray}(e1, e2)
\S4method{pmax2}{SC_GDSArray,SC_GDSArray}(e1, e2)
\S4method{pmax2}{SC_GDSArray,vector}(e1, e2)
\S4method{pmax2}{vector,SC_GDSArray}(e1, e2)
}

\arguments{
    \item{x}{A \link{SC_GDSArray} or \link{SC_GDSMatrix} object}
    \item{i, j, ...}{indices specifying elements to extract}
    \item{drop}{if \code{TRUE} the result will be coerced to the lowest
        possible dimension}
    \item{e1, e2}{objects}
    \item{value}{\code{NULL}, a character vector for \code{names<-} or a list
        of character vectors for \code{dimnames<-}}
    \item{center}{either a logical value or a numeric vector (e.g., \code{FALSE}
        or \code{0} for no centering)}
    \item{scale}{either a logical value or a numeric vector (e.g., \code{TRUE}
        or \code{1} for no scaling)}
}
\value{
    All these operations return a SC_GDSArray or SC_GDSMatrix object.
}

\details{
    All these operations return a SC_GDSArray or SC_GDSMatrix object.
    \describe{
        \item{\code{Arith}:}{
            \code{"+"}, \code{"-"}, \code{"*"}, \code{"^"}, \code{"\%\%"},
            \code{"\%/\%"}, \code{"/"}}
        \item{\code{Compare}:}{
            \code{"=="}, \code{">"}, \code{"<"}, \code{"!="}, \code{"<="},
            \code{">="}}
        \item{\code{Logic}:}{
            \code{"&"}, \code{"|"}.}
        \item{\code{Ops}:}{
            \code{"Arith"}, \code{"Compare"}, \code{"Logic"}}
        \item{\code{Math}:}{
            \code{"abs"}, \code{"sign"}, \code{"sqrt"}, \code{"ceiling"},
            \code{"floor"}, \code{"trunc"}, \code{"cummax"}, \code{"cummin"},
            \code{"cumprod"}, \code{"cumsum"}, \code{"log"}, \code{"log10"},
            \code{"log2"}, \code{"log1p"}, \code{"acos"}, \code{"acosh"},
            \code{"asin"}, \code{"asinh"}, \code{"atan"}, \code{"atanh"},
            \code{"exp"}, \code{"expm1"},
            \code{"cos"}, \code{"cosh"}, \code{"cospi"},
            \code{"sin"}, \code{"sinh"}, \code{"sinpi"},
            \code{"tan"}, \code{"tanh"}, \code{"tanpi"},
            \code{"gamma"}, \code{"lgamma"}, \code{"digamma"},
            \code{"trigamma"}
        }
    }
}

\author{Xiuwen Zheng}
\seealso{
    \link{Ops}, \link{Math}, \link{SCArray-stats}
}

\examples{
fn <- system.file("extdata", "example.gds", package="SCArray")

x <- scArray(fn, "counts")

x[1:8, 1:32]
x > 0
pmin2(x, 1)
log1p(x)
scale(x)

rm(x)
}

\keyword{methods}
\keyword{GDS}
