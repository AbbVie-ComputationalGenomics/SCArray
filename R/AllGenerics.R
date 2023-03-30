#######################################################################
#
# Package name: SCArray
#
# Description:
#     Large-scale single-cell RNA-seq data manipulation with GDS files
#
# Copyright (C) 2021-2023    Xiuwen Zheng (@AbbVie-ComputationalGenomics)
# License: GPL-3
#


#######################################################################

DClass  <- "DelayedArray"
DMatrix <- "DelayedMatrix"
SClass  <- "SC_GDSArray"
SMatrix <- "SC_GDSMatrix"


# For internal use only
setMethod("matrixClass", SClass, function(x) SMatrix)

# Automatic coercion method from GDSArray to GDSMatrix (muted for
# higher dimensions) this function works only when GDSArray is
# 2-dimensional, otherwise it fails.

setAs(SClass, SMatrix, function(from) new(SMatrix, from))
setAs(SMatrix, SClass, function(from) from)
setAs("ANY", SMatrix,
    function(from) as(as(from, SClass), SMatrix))

# setMethod(
#     "DelayedArray", "SCArraySeed",
#     function(seed) new_DelayedArray(seed, Class=SClass)
# )


#######################################################################

.reset_OPS_env <- function(x, varnm)
{
    # find OPS
    if (is(x, "DelayedUnaryIsoOpStack"))
    {
        for (f in x@OPS)
        {
            e <- environment(f)
            if (!is.null(e))
            {
                for (nm in varnm)
                    if (exists(nm, envir=e)) remove(list=nm, envir=e)
            }
        }
        return(TRUE)
    }
    if (is(x, "DelayedUnaryOp"))
        return(.reset_OPS_env(x@seed, varnm))
    # find OP
    if (is(x, "DelayedNaryOp"))
    {
        if (is(x, "DelayedNaryIsoOp"))
        {
            e <- environment(x@OP)
            if (!is.null(e))
            {
                for (nm in varnm)
                    if (exists(nm, envir=e)) remove(list=nm, envir=e)
            }
            return(TRUE)
        }
        x <- x@seeds
    }
    # find next
    if (is.list(x) && !is.array(x))
    {
        for (y in x)
        {
            if (.reset_OPS_env(y, varnm)) return(TRUE)
        }
    }
    FALSE
}


.sc_val <- function(v)
{
    if (is(v, DMatrix) && !is(v, SMatrix))
        as(v, SMatrix)
    else if (is(v, DClass) && !is(v, SClass))
        as(v, SClass)
    else
        v
}

.sc_val_e1 <- function(v)
{
    .reset_OPS_env(v, "e1")
    .sc_val(v)
}

.sc_val_e2 <- function(v)
{
    .reset_OPS_env(v, "e2")
    .sc_val(v)
}



# subsetting
setMethod("[", SClass,
    function(x, i, j, ..., drop=TRUE) .sc_val(callNextMethod()) )
setMethod("[[", SClass,
    function(x, i, j, ...) .sc_val(callNextMethod()) )

# transpose
setMethod("aperm", SClass,
    function(a, perm, ...) .sc_val(callNextMethod()) )

# names<-
setMethod("names<-", SClass, function(x, value) .sc_val(callNextMethod()) )

# dimnames<-
setMethod("dimnames<-", SClass, function(x, value) .sc_val(callNextMethod()) )


# Ops
setMethod("Ops", c(SClass, "vector"),
    function(e1, e2) .sc_val_e1(callGeneric(as(e1, DClass), e2)) )
setMethod("Ops", c("vector", SClass),
    function(e1, e2) .sc_val_e2(callGeneric(e1, as(e2, DClass))) )
setMethod("Ops", c(SClass, SClass),
    function(e1, e2) .sc_val(callGeneric(as(e1, DClass), as(e2, DClass))) )
setMethod("+", c(SClass, "missing"),    # unary operators "+"
    function(e1, e2) .sc_val_e1(callGeneric(as(e1, DClass))) )
setMethod("-", c(SClass, "missing"),    # unary operators "-"
    function(e1, e2) .sc_val_e1(callGeneric(as(e1, DClass))) )


# Math
setMethods("Math", SClass,
    function(x) .sc_val(callGeneric(as(x, DClass))) )

# pmax2
setMethod("pmax2", c(SClass, "vector"),
    function(e1, e2) .sc_val_e1(callGeneric(as(e1, DClass), e2)) )
setMethod("pmax2", c("vector", SClass),
    function(e1, e2) .sc_val_e2(callGeneric(e1, as(e2, DClass))) )
setMethod("pmax2", c(SClass, SClass),
    function(e1, e2) .sc_val(callGeneric(as(e1, DClass), as(e2, DClass))) )

# pmin2
setMethod("pmin2", c(SClass, "vector"),
    function(e1, e2) .sc_val_e1(callGeneric(as(e1, DClass), e2)) )
setMethod("pmin2", c("vector", SClass),
    function(e1, e2) .sc_val_e2(callGeneric(e1, as(e2, DClass))) )
setMethod("pmin2", c(SClass, SClass),
    function(e1, e2) .sc_val(callGeneric(as(e1, DClass), as(e2, DClass))) )


# scale for each column
x_scale <- function(x, center=TRUE, scale=TRUE)
{
    # check
    x_check(x, "Calling SCArray:::x_scale() with %s ...")
    stopifnot(is(x, "SC_GDSMatrix"))
    stopifnot(is.logical(center) || is.numeric(center),
        length(center)==1L || length(center)==ncol(x))
    stopifnot(is.logical(scale) || is.numeric(scale),
        length(scale)==1L || length(scale)==ncol(x))
    # get center and/or scale
    if (isTRUE(center) && isTRUE(scale))
    {
        v <- scColMeanVar(x, na.rm=TRUE)  # calculate mean and var together
        center <- v[,1L]
        scale <- sqrt(v[,2L])
    } else if (isTRUE(center))
    {
        center <- colMeans(x)
    } else if (isTRUE(scale))
    {
        scale <- rowSds(x)
    }
    center <- rep(as.numeric(center), length.out=ncol(x))
    if (isFALSE(scale)) scale <- 1
    scale <- rep(as.numeric(scale), length.out=ncol(x))
    # do
    x <- t(x)
    attr_c <- !isTRUE(all(center==0))
    attr_s <- !isTRUE(all(scale==1))
    if (attr_c) x <- x - center
    if (attr_s)
    {
        inv <- 1/scale
        inv[!is.finite(inv)] <- 0
        x <- x * inv
    }
    x <- t(x)
    if (attr_c) attr(x, "scaled:center") <- center
    if (attr_s) attr(x, "scaled:scale") <- scale
    # output
    x
}

setMethod("scale", SMatrix, x_scale)


#######################################################################
# Generic functions -- scGetFiles()

setGeneric("scGetFiles", function(object, ...) standardGeneric("scGetFiles"))

setMethod("scGetFiles", "SC_GDSArray", function(object, ...) {
    s <- seedApply(object, function(x)
        tryCatch(path(x), error=function(e) NULL)
    )
    unique(unlist(s))
})

setMethod("scGetFiles", "SummarizedExperiment",
    function(object, ...) unique(unlist(lapply(assays(object), scGetFiles))))


#######################################################################
# Generic functions -- scMemory()

setGeneric("scMemory", function(x, ...) x)

setMethod("scMemory", "DelayedArray", function(x, ...)
{
    if (is_sparse(x) && length(dim(x))==2L)
        as(x, "sparseMatrix")
    else
        as.array(x)
})

setMethod("scMemory", "SummarizedExperiment", function(x, ...)
{
    for (i in seq_along(assays(x)))
        assays(x)[[i]] <- scMemory(assays(x)[[i]])
    x
})


#######################################################################
# Generic functions -- scMemory()

scRowAutoGrid <- function(x, force=FALSE, nnzero=NULL)
{
    # check
    stopifnot(is.matrix(x) || is(x, "Matrix") || is(x, "DelayedMatrix"))
    stopifnot(is.logical(force))
    if (!is.null(nnzero))
    {
        stopifnot(is.numeric(nnzero), length(nnzero)==nrow(x))
        if (is.double(nnzero)) nnzero <- as.integer(nnzero)
        z <- nnzero < 0L
        if (anyNA(z)) stop("'nnzero' should have no NA.")
        if (any(z)) stop("'nnzero' should be non-negative.")
    }
    # do
    ans <- rowAutoGrid(x)
    as.sparse <- NA
    is_sp <- is(x, "sparseMatrix")
    if ((is(x, SMatrix) && is_sparse(x)) || is_sp)
    {
        if (isTRUE(force) || is_sp || x_type(x)==1L)
        {
            # x_type(x)==1L (not transposed, column-oriented)
            # the number of non-zeros for each row
            if (is.null(nnzero))
                nnzero <- as.integer(row_nnzero(x, na.counted=TRUE))
            # to use SafeArrayViewport
            bs <- .Call(c_sparse_blocksize, getAutoBlockSize(),
                .Machine$integer.max / ncol(x), nnzero, integer(length(nnzero)))
            v <- ArbitraryArrayGrid(list(cumsum(bs), ncol(x)))
            as.sparse <- isTRUE(force) || length(v)<=length(ans)
            if (as.sparse) ans <- v
        }
    }
    # output
    attr(ans, "as.sparse") <- as.sparse
    ans
}

scColAutoGrid <- function(x, force=FALSE, nnzero=NULL)
{
    # check
    stopifnot(is.matrix(x) || is(x, "Matrix") || is(x, "DelayedMatrix"))
    stopifnot(is.logical(force))
    if (!is.null(nnzero))
    {
        stopifnot(is.numeric(nnzero), length(nnzero)==ncol(x))
        if (is.double(nnzero)) nnzero <- as.integer(nnzero)
        z <- nnzero < 0L
        if (anyNA(z)) stop("'nnzero' should have no NA.")
        if (any(z)) stop("'nnzero' should be non-negative.")
    }
    # do
    ans <- colAutoGrid(x)
    as.sparse <- NA
    is_sp <- is(x, "sparseMatrix")
    if ((is(x, SMatrix) && is_sparse(x)) || is_sp)
    {
        if (isTRUE(force) || is_sp || x_type(x)==2L)
        {
            # x_type(x)==2L (transposed, row-oriented)
            # the number of non-zeros for each column
            if (is.null(nnzero))
                nnzero <- as.integer(col_nnzero(x, na.counted=TRUE))
            # to use SafeArrayViewport
            bs <- .Call(c_sparse_blocksize, getAutoBlockSize(),
                .Machine$integer.max / nrow(x), nnzero, integer(length(nnzero)))
            v <- ArbitraryArrayGrid(list(nrow(x), cumsum(bs)))
            as.sparse <- isTRUE(force) || length(v)<=length(ans)
            if (as.sparse) ans <- v
        }
    }
    # output
    attr(ans, "as.sparse") <- as.sparse
    ans
}


#######################################################################
# Internal functions

.get_num_worker <- function(BPPARAM)
{
    if (is.null(BPPARAM))
        1L
    else
        bpnworkers(BPPARAM)
}

####

.parallel_col_apply <- function(x, BPPARAM, Fun, as.sparse=NA, .flatten=TRUE,
    .progress=NA, ...)
{
    stopifnot(is.null(BPPARAM) || is(BPPARAM, "BiocParallelParam"))
    # progress bar
    stopifnot(is.logical(.progress), length(.progress)==1L)
    if (is.na(.progress)) .progress <- x_progress_verbose()
    # split columns
    sp <- scNumSplit(ncol(x), BPPARAM)
    if (length(sp) > 1L)
    {
        x_msg(sprintf("\\=> Distributed to %d processes ...", length(sp)))
        # distribute
        if (!is.null(dimnames(x))) dimnames(x) <- NULL  # no need dimnames
        lst <- bplapply(sp, function(s, Fun, x, as.sparse, ...)
        {
            # sub columns
            y <- x[, seq.int(s[1L], s[2L]), drop=FALSE]
            # reduce the sub matrix
            blockApply(y, Fun, grid=colAutoGrid(y), BPPARAM=NULL,
                as.sparse=as.sparse, ...)
        }, BPPARAM=BPPARAM, Fun=Fun, x=x, as.sparse=as.sparse, ...)
        # reduce
        ans <- unlist(lst, recursive=FALSE, use.names=FALSE)
    } else {
        # sequentially
        gd <- colAutoGrid(x)
        if (.progress)
        {
            pb <- txtProgressBar(0L, length(gd), style=3L, width=64L,
                file=stderr())
            on.exit(close(pb))
            ans <- blockApply(x, function(b, .fun, .pb, ...)
                {
                    setTxtProgressBar(.pb, currentBlockId())
                    .fun(b, ...)
                }, grid=gd, as.sparse=as.sparse, BPPARAM=NULL,
                    .fun=Fun, .pb=pb, ...)
        } else {
            ans <- blockApply(x, Fun, grid=gd, as.sparse=as.sparse,
                BPPARAM=NULL, ...)
        }
    }
    # output
    if (.flatten) ans <- unlist(ans, recursive=TRUE, use.names=FALSE)
    ans
}

####

.parallel_col_reduce <- function(x, BPPARAM, Fun, InitFun, ReduceFun,
    as.sparse=NA, .progress=NA, ...)
{
    stopifnot(is.null(BPPARAM) || is(BPPARAM, "BiocParallelParam"))
    # progress bar
    stopifnot(is.logical(.progress), length(.progress)==1L)
    if (is.na(.progress)) .progress <- x_progress_verbose()
    # split columns
    sp <- scNumSplit(ncol(x), BPPARAM)
    if (length(sp) > 1L)
    {
        x_msg(sprintf("\\=> Distributed to %d processes ...", length(sp)))
        # distribute
        if (!is.null(dimnames(x))) dimnames(x) <- NULL  # no need dimnames
        lst <- bplapply(sp, function(s, Fun, x, as.sparse, ...)
        {
            # sub columns
            y <- x[, seq.int(s[1L], s[2L]), drop=FALSE]
            attr(y, "col_range") <- s
            # initial value
            if (is.function(InitFun))
                init <- InitFun(y, ...)
            else
                init <- InitFun
            # reduce the sub matrix
            blockReduce(Fun, y, init, grid=colAutoGrid(y), as.sparse=as.sparse,
                ...)
        }, BPPARAM=BPPARAM, Fun=Fun, x=x, as.sparse=as.sparse, ...)
        # reduce
        base::Reduce(ReduceFun, lst)
    } else {
        # sequentially
        if (is.function(InitFun))
            init <- InitFun(x, ...)
        else
            init <- InitFun
        gd <- colAutoGrid(x)
        if (.progress)
        {
            pb <- txtProgressBar(0L, length(gd), style=3L, width=64L,
                file=stderr())
            on.exit(close(pb))
            blockReduce(function(b, v, .fun, .pb, ...)
                {
                    setTxtProgressBar(.pb, currentBlockId())
                    .fun(b, v, ...)
                }, x, init, grid=gd, as.sparse=as.sparse, .fun=Fun, .pb=pb, ...)
        } else {
            blockReduce(Fun, x, init, grid=gd, as.sparse=as.sparse, ...)
        }
    }
}

####

.parallel_col_reduce2 <- function(x, BPPARAM, Fun, InitFun, ReduceFun,
    as.sparse=NA, .progress=NA, ...)
{
    stopifnot(is.null(BPPARAM) || is(BPPARAM, "BiocParallelParam"))
    # progress bar
    stopifnot(is.logical(.progress), length(.progress)==1L)
    if (is.na(.progress)) .progress <- x_progress_verbose()
    # split columns
    sp <- scNumSplit(ncol(x), BPPARAM)
    if (length(sp) > 1L)
    {
        x_msg(sprintf("\\=> Distributed to %d processes ...", length(sp)))
        # distribute
        if (!is.null(dimnames(x))) dimnames(x) <- NULL  # no need dimnames
        lst <- bplapply(sp, function(s, Fun, x, as.sparse, ...)
        {
            # sub columns
            y <- x[, seq.int(s[1L], s[2L]), drop=FALSE]
            attr(y, "col_range") <- s
            # initial value
            if (is.function(InitFun))
                init <- InitFun(y, split=s, ...)
            else
                init <- InitFun
            # reduce the sub matrix
            blockReduce(Fun, y, init, grid=colAutoGrid(y), as.sparse=as.sparse,
                split=s, ...)
        }, BPPARAM=BPPARAM, Fun=Fun, x=x, as.sparse=as.sparse, ...)
        # reduce
        base::Reduce(ReduceFun, lst)
    } else {
        # sequentially
        sp <- c(1L, ncol(x))
        if (is.function(InitFun))
            init <- InitFun(x, split=sp, ...)
        else
            init <- InitFun
        gd <- colAutoGrid(x)
        if (.progress)
        {
            pb <- txtProgressBar(0L, length(gd), style=3L, width=64L,
                file=stderr())
            on.exit(close(pb))
            blockReduce(function(b, v, .fun, .pb, ...)
            {
                setTxtProgressBar(.pb, currentBlockId())
                .fun(b, v, ...)
            }, x, init, grid=gd, as.sparse=as.sparse, .fun=Fun, .pb=pb,
                split=sp, ...)
        } else {
            blockReduce(Fun, x, init, grid=gd, as.sparse=as.sparse, split=sp,
                ...)
        }
    }
}


################

x_verbose <- function()
{
    isTRUE(getOption("SCArray.verbose"))
}

x_progress_verbose <- function()
{
    isTRUE(getOption("SCArray.verbose")) ||
    isTRUE(getOption("SCArray.progress.verbose"))
}

x_progress <- function(grid, verbose=TRUE)
{
    if (verbose && x_progress_verbose())
        txtProgressBar(0L, length(grid), style=3L, width=64L, file=stderr())
    else
        NULL
}

x_msg <- function(msg)
{
    if (x_verbose()) message(msg)
}

x_check <- function(x, msg)
{
    if (!is.null(x))
    {
        if (x_verbose())
        {
            if (grepl("%s", msg, fixed=TRUE))
            {
                s <- class(x)[1L]
                if (x_type(x) == 2L) s <- paste("transposed", s)
                s <- paste0(s, " [", paste(dim(x), collapse="x"), "]")
                msg <- sprintf(msg, s)
            }
            message(msg)
        }
        DelayedArray:::.get_ans_type(x, must.be.numeric=TRUE)
        stopifnot(length(dim(x)) == 2L)
    } else {
        if (x_verbose()) message(msg)
    }
    invisible()
}

x_subset <- function(x, rows, cols)
{
    if (!is.null(rows) && !is.null(cols))
    {
        x <- x[rows, cols, drop = FALSE]
    } else if (!is.null(rows))
    {
        x <- x[rows, , drop = FALSE]
    } else if (!is.null(cols))
    {
        x <- x[, cols, drop = FALSE]
    }
    x
}

x_nperm <- function(x)
{
    ans <- if (is(x, "DelayedAperm")) 1L else 0L
    if (is(x, "DelayedUnaryOp"))
        return(ans + x_nperm(x@seed))
    if (is(x, "DelayedNaryOp"))
        x <- x@seeds
    if (is.list(x) && !is.array(x))
    {
        for (y in x) ans <- ans + x_nperm(y)
    }
    ans
}

# Return 1 for SCArraySeed, 2L for transposed SCArraySeed, 3L for others
x_type <- function(x)
{
    if (nseed(x)==1L && is(seed(x), "SCArraySeed"))
    {
        # whether transposed or not
        if (x_nperm(x) %% 2L == 0L) 1L else 2L
    } else
        3L
}


