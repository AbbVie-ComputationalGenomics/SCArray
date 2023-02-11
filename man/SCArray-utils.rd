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
}

\arguments{
    \item{x}{A \link{SC_GDSArray} or \link{SC_GDSMatrix} object}
    \item{i, j, ...}{indices specifying elements to extract}
    \item{drop}{if \code{TRUE} the result will be coerced to the lowest
        possible dimension}
    \item{e1, e2}{objects}
    \item{value}{\code{NULL}, a character vector for \code{names<-} or a list
        of character vectors for \code{dimnames<-}}
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
1
}

\keyword{methods}
\keyword{GDS}
