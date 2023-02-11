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
    if (is_sparse(x))
        as(x, "sparseMatrix")
    else
        as.matrix(x)
})

setMethod("scMemory", "SummarizedExperiment", function(x, ...)
{
    for (i in seq_along(assays(x)))
        assays(x)[[i]] <- scMemory(assays(x)[[i]])
    x
})



#######################################################################
# Internal functions

x_verbose <- function()
{
    isTRUE(getOption("SCArray.verbose"))
}

x_progress_verbose <- function()
{
    isTRUE(getOption("SCArray.verbose")) ||
    isTRUE(getOption("SCArray.progress.verbose"))
}

x_progress <- function(grid)
{
    if (x_progress_verbose())
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


