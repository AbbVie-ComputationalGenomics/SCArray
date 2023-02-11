#######################################################################
#
# Package name: SCArray
#
# Description:
#     Large-scale single-cell RNA-seq data manipulation with GDS files
#
# Copyright (C) 2020-2023    Xiuwen Zheng (@AbbVie-ComputationalGenomics)
# License: GPL-3
#


#######################################################################
# Class definition

setClass("SCArrayFileClass", contains="gds.class")

# test the validity of objects
setValidity("SCArrayFileClass", function(object)
    {
        if (!inherits(object, "gds.class"))
            return("'object' should be inherited from 'gds.class'.")
        TRUE
    }
)

# to avoid dropping the names in S3 class 'gds.class'
setMethod("updateObject", "SCArrayFileClass",
    function(object, ..., verbose=FALSE) object)



#######################################################################
# Class definition
# See the vignette of the DelayedArray package
#

# if the file is open, no action internally
.reopen <- function(x) gdsfmt:::.reopen(x@gds)

# slot member access
.gds      <- function(x) x@gds
.filename <- function(x) x@filename
.varname  <- function(x) x@varname
.dim      <- function(x) x@dim


# define SCArraySeed class
setClass("SCArraySeed", contains="Array",
    slots = c(
        gds = "SCArrayFileClass",
        filename = "character",
        varname = "character",
        dim = "integer"
    )
)

setClass("SC_GDSArray", contains="DelayedArray")

setClass("SC_GDSMatrix", contains=c("DelayedMatrix", "SC_GDSArray"))



# set the DelayedArray function
setMethod("DelayedArray", "SCArraySeed",
    function(seed) new_DelayedArray(seed, Class="SC_GDSArray")
)


# show method for SCArraySeed object
setMethod("show", "SCArraySeed", function(object)
    {
        .cat("SCArraySeed\n",
            "File: ", .filename(object), "\n",
            "Array node: ", .varname(object), "\n",
            "Dim: ", paste(.dim(object), collapse=" x "))
    }
)


# extract a dense array from DelayedArray
.extract_array_sc_seed <- function(x, index)
{
    # check
    stopifnot(is.list(index), !anyNA(index))
    ans_dim <- DelayedArray:::get_Nindex_lengths(index, dim(x))
    # reopen the file if needed
    .reopen(x)
    # read
    if (any(ans_dim == 0L))
    {
        tp <- objdesp.gdsn(index.gdsn(.gds(x), .varname(x)))$type
        ans <- switch(as.character(tp),
            Raw=raw(), Integer=integer(), Logical=logical(),
            Real=double(), String=character(),
            stop("Unsupported data type: ", tp))
        dim(ans) <- ans_dim
    } else {
        nd <- index.gdsn(.gds(x), .varname(x))
        ans <- readex.gdsn(nd, index, .sparse=FALSE)
        if (!is.array(ans))  # ans must be an array
            dim(ans) <- ans_dim
    }
    ans
}

setMethod("extract_array", "SCArraySeed", .extract_array_sc_seed)


# return whether the array is sparse or not
.is_sparse_sc_seed <- function(x)
{
    # reopen the file if needed
    .reopen(x)
    # whether is a sparse array
    is.sparse.gdsn(index.gdsn(x@gds, x@varname))
}

setMethod("is_sparse", "SCArraySeed", .is_sparse_sc_seed)


# extract a sparse array from DelayedArray
.extract_sparse_sc_seed <- function(x, index)
{
    # check
    stopifnot(is.list(index), !anyNA(index))
    ans_dim <- DelayedArray:::get_Nindex_lengths(index, dim(x))
    # reopen the file if needed
    .reopen(x)
    # read
    if (any(ans_dim == 0L))
    {
        tp <- objdesp.gdsn(index.gdsn(.gds(x), .varname(x)))$type
        ans <- switch(as.character(tp),
            Raw=raw(), Integer=integer(), Logical=logical(),
            Real=double(), String=character(),
            stop("Unsupported data type: ", tp))
        SparseArraySeed(ans_dim, nzdata=ans)
    } else {
        nd <- index.gdsn(.gds(x), .varname(x))
        ans <- readex.gdsn(nd, index, .sparse=TRUE)
        if (inherits(ans, "gds_sparse_nz_class"))
        {
            # ans is a list(nzdata, nzindex)
            SparseArraySeed(ans_dim, ans$nzindex, ans$nzdata, check=FALSE)
        } else {
            # ans could be a dense array, dgCMatrix or lgCMatrix
            as(ans, "SparseArraySeed")
        }
    }
}

setMethod("extract_sparse_array", "SCArraySeed", .extract_sparse_sc_seed)


# SCArraySeed constructor
SCArraySeed <- function(gds, varname)
{
    # check gds
    stopifnot(inherits(gds, "SCArrayFileClass"))
    # check varname
    stopifnot(is.character(varname), length(varname)==1L, !is.na(varname))
    nd <- index.gdsn(gds, varname, silent=TRUE)
    if (is.null(nd))
        stop("No '", varname, "'")
    # check dimension
    dp <- objdesp.gdsn(nd)
    if (!dp$is.array)
        stop("'", varname, "' is not an array.")
    dm <- dp$dim
    # output
    new2("SCArraySeed", gds=gds, filename=gds$filename, varname=varname, dim=dm)
}


setMethod("path", "SCArraySeed", function(object) object@filename)
setMethod("path", "SCArrayFileClass", function(object) object$filename)
