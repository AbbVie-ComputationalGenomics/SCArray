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


################

.x_row_sums <- function(x, na.rm, ...)
{
    blockReduce(function(bk, v, na.rm) {
        .Call(c_rowSums, bk, v, na.rm)
    }, x, double(nrow(x)), grid=colAutoGrid(x), as.sparse=NA, na.rm=na.rm, ...)
}

.x_col_sums <- function(x, na.rm, ...)
{
    unlist(blockApply(x, function(bk, na.rm) {
        .Call(c_colSums, bk, na.rm)
    }, grid=colAutoGrid(x), as.sparse=NA, na.rm=na.rm, ...))
}

x_rowSums <- function(x, na.rm=FALSE, dims=1)
{
    x_check(x, "Calling SCArray:::x_rowSums() with %s ...")
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    stopifnot(identical(dims, 1))
    # output
    switch(x_type(x),
        .x_row_sums(x, na.rm),     # 1
        .x_col_sums(t(x), na.rm),  # 2
        {   # 3
            x_msg("Calling DelayedArray::rowSums() ...")
            callNextMethod()
        }
    )
}

x_colSums <- function(x, na.rm=FALSE, dims=1)
{
    x_check(x, "Calling SCArray:::x_colSums() with %s ...")
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    stopifnot(identical(dims, 1))
    # output
    switch(x_type(x),
        .x_col_sums(x, na.rm),     # 1
        .x_row_sums(t(x), na.rm),  # 2
        {   # 3
            x_msg("Calling DelayedArray::colSums() ...")
            callNextMethod()
        }
    )
}

x_rowSums2 <- function(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
{
    x_check(x, "Calling SCArray:::x_rowSums2() with %s ...")
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    k <- x_type(x)
    if (k < 3L)
    {
        x <- x_subset(x, rows, cols)
        if (k  == 1L)
            v <- .x_row_sums(x, na.rm, ...)
        else
            v <- .x_col_sums(t(x), na.rm, ...)
        if (isTRUE(useNames)) names(v) <- rownames(x)
        v  # output
    } else {
        x_msg("Calling DelayedMatrixStats::rowSums2() ...")
        callNextMethod()
    }
}

x_colSums2 <- function(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
{
    x_check(x, "Calling SCArray:::x_colSums2() with %s ...")
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    k <- x_type(x)
    if (k < 3L)
    {
        x <- x_subset(x, rows, cols)
        if (k  == 1L)
            v <- .x_col_sums(x, na.rm, ...)
        else
            v <- .x_row_sums(t(x), na.rm, ...)
        if (isTRUE(useNames)) names(v) <- rownames(x)
        v  # output
    } else {
        x_msg("Calling DelayedMatrixStats::colSums2() ...")
        callNextMethod()
    }
}

setMethod("rowSums", SMatrix, x_rowSums)
setMethod("colSums", SMatrix, x_colSums)
setMethod("rowSums2", SMatrix, x_rowSums2)
setMethod("colSums2", SMatrix, x_colSums2)


################

.x_row_prods <- function(x, na.rm, ...)
{
    blockReduce(function(bk, v, na.rm) {
        .Call(c_rowProds, bk, v, na.rm)
    }, x, init=rep(1, nrow(x)), grid=colAutoGrid(x), na.rm=na.rm, ...)
}

.x_col_prods <- function(x, na.rm, ...)
{
    unlist(blockApply(x, function(bk, na.rm) {
        .Call(c_colProds, bk, na.rm)
    }, grid=colAutoGrid(x), na.rm=na.rm, ...))
}

x_rowProds <- function(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
{
    x_check(x, "Calling SCArray:::x_rowProds() with %s ...")
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    if (x_type(x) == 1L)
    {
        x <- x_subset(x, rows, cols)
        v <- .x_row_prods(x, na.rm, ...)
        if (isTRUE(useNames)) names(v) <- rownames(x)
        v  # output
    } else {
        x_msg("Calling DelayedMatrixStats::rowProds() ...")
        callNextMethod()
    }
}

x_colProds <- function(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
{
    x_check(x, "Calling SCArray:::x_colProds() with %s ...")
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    if (x_type(x) != 2L)
    {
        x_msg("Calling DelayedMatrixStats::colProds() ...")
        callNextMethod()
    } else {
        x_msg("Calling DelayedArray::.x_row_prods() ...")
        x <- x_subset(x, rows, cols)
        v <- .x_row_prods(t(x), na.rm, ...)
        if (isTRUE(useNames)) names(v) <- colnames(x)
        v  # output
    }
}

setMethod("rowProds", SMatrix, x_rowProds)
setMethod("colProds", SMatrix, x_colProds)


################

.x_row_means <- function(x, na.rm, ...)
{
    # block read
    rv <- blockReduce(function(bk, v, na.rm) {
        .Call(c_rowMeans, bk, v, na.rm)
    }, x, double(nrow(x)*2L), grid=colAutoGrid(x), as.sparse=NA, na.rm=na.rm, ...)
    # finally
    .Call(c_rowMeans_final, rv)
}

.x_col_means <- function(x, na.rm, ...)
{
    unlist(blockApply(x, function(bk, na.rm) {
        .Call(c_colMeans, bk, na.rm)
    }, grid=colAutoGrid(x), as.sparse=NA, na.rm=na.rm, ...))
}

x_rowMeans <- function(x, na.rm=FALSE, dims=1)
{
    x_check(x, "Calling SCArray:::x_rowMeans() with %s ...")
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    stopifnot(identical(dims, 1))
    # output
    switch(x_type(x),
        .x_row_means(x, na.rm),     # 1
        .x_col_means(t(x), na.rm),  # 2
        {   # 3
            x_msg("Calling DelayedArray::rowMeans() ...")
            callNextMethod()
        }
    )
}

x_colMeans <- function(x, na.rm=FALSE, dims=1)
{
    x_check(x, "Calling SCArray:::x_colMeans() with %s ...")
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    stopifnot(identical(dims, 1))
    switch(x_type(x),
        .x_col_means(x, na.rm),     # 1
        .x_row_means(t(x), na.rm),  # 2
        {   # 3
            x_msg("Calling DelayedArray::colMeans() ...")
            callNextMethod()
        }
    )
}

x_rowMeans2 <- function(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
{
    x_check(x, "Calling SCArray:::x_rowMeans2() with %s ...")
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    k <- x_type(x)
    if (k < 3L)
    {
        x <- x_subset(x, rows, cols)
        if (k == 1L)
            v <- .x_row_means(x, na.rm, ...)
        else
            v <- .x_col_means(t(x), na.rm, ...)
        if (isTRUE(useNames)) names(v) <- rownames(x)
        v
    } else {
        x_msg("Calling DelayedMatrixStats::rowMeans2() ...")
        callNextMethod()
    }
}

x_colMeans2 <- function(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
{
    x_check(x, "Calling SCArray:::x_colMeans2() with %s ...")
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    k <- x_type(x)
    if (k < 3L)
    {
        x <- x_subset(x, rows, cols)
        if (k == 1L)
            v <- .x_col_means(x, na.rm, ...)
        else
            v <- .x_row_means(t(x), na.rm, ...)
        if (isTRUE(useNames)) names(v) <- colnames(x)
        v
    } else {
        x_msg("Calling DelayedMatrixStats::colMeans2() ...")
        callNextMethod()
    }
}

setMethod("rowMeans", SMatrix, x_rowMeans)
setMethod("colMeans", SMatrix, x_colMeans)
setMethod("rowMeans2", SMatrix, x_rowMeans2)
setMethod("colMeans2", SMatrix, x_colMeans2)


################

.x_row_w_means <- function(x, w, na.rm, ...)
{
    # initialize
    stopifnot(is.numeric(w), length(w)==ncol(x))
    if (is.integer(w)) w <- as.double(w)
    .Call(c_init_block)
    # block read
    rv <- blockReduce(function(bk, v, w, na.rm) {
        .Call(c_rowWMeans, bk, v, w, na.rm)
    }, x, double(nrow(x)*2L), grid=colAutoGrid(x), as.sparse=NA, w=w, na.rm=na.rm)
    # finally
    .Call(c_rowWMeans_final, rv)
}

.x_col_w_means <- function(x, w, na.rm, ...)
{
    stopifnot(is.numeric(w), length(w)==nrow(x))
    if (is.integer(w)) w <- as.double(w)
    # block read
    unlist(blockApply(x, function(bk, w, na.rm) {
        .Call(c_colMeans, bk, w, na.rm)
    }, grid=colAutoGrid(x), as.sparse=NA, w=w, na.rm=na.rm, ...))
}

x_rowWeightedMeans <- function(x, w=NULL, rows=NULL, cols=NULL, na.rm=FALSE,
    ..., useNames=NA)
{
    x_check(x, "Calling SCArray:::x_rowWeightedMeans() with %s ...")
    if (is.null(w))
        return(x_rowMeans2(x, rows, cols, na.rm, ..., useNames=useNames))
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    k <- x_type(x)
    if (k < 3L)
    {
        x <- x_subset(x, rows, cols)
        if (k == 1L)
            v <- .x_row_w_means(x, w, na.rm, ...)
        else
            v <- .x_col_w_means(t(x), w, na.rm, ...)
        if (isTRUE(useNames)) names(v) <- rownames(x)
        v
    } else {
        x_msg("Calling DelayedMatrixStats::rowWeightedMeans() ...")
        callNextMethod()
    }
}

x_colWeightedMeans <- function(x, w=NULL, rows=NULL, cols=NULL, na.rm=FALSE,
    ..., useNames=NA)
{
    x_check(x, "Calling SCArray:::x_colWeightedMeans() with %s ...")
    if (is.null(w))
        return(x_colMeans2(x, rows, cols, na.rm, ..., useNames=useNames))
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    k <- x_type(x)
    if (k < 3L)
    {
        x <- x_subset(x, rows, cols)
        if (k == 1L)
            v <- .x_col_w_means(x, w, na.rm, ...)
        else
            v <- .x_row_w_means(t(x), w, na.rm, ...)
        if (isTRUE(useNames)) names(v) <- colnames(x)
        v
    } else {
        x_msg("Calling DelayedMatrixStats::colWeightedMeans() ...")
        callNextMethod()
    }
}

setMethod("rowWeightedMeans", SMatrix, x_rowWeightedMeans)
setMethod("colWeightedMeans", SMatrix, x_colWeightedMeans)


################

.x_row_vars <- function(x, na.rm, center, ...)
{
    # check
    if (length(center))
    {
        if (length(center)==1L) center <- rep(center, nrow(x))
        stopifnot(nrow(x) == length(center))
        if (is.integer(center)) center <- as.double(center)
    }
    # block read
    v <- blockReduce(function(bk, v, na.rm, center) {
        .Call(c_rowVars, bk, v, na.rm, center)
    }, x, init=double(nrow(x)*3L), grid=colAutoGrid(x), as.sparse=NA,
        na.rm=na.rm, center=center, ...)
    # finally
    .Call(c_rowVars_final, v, center)
}

.x_col_vars <- function(x, na.rm, center, ...)
{
    # check
    if (length(center))
    {
        if (length(center)==1L) center <- rep(center, ncol(x))
        stopifnot(ncol(x) == length(center))
        if (is.integer(center)) center <- as.double(center)
    }
    # block read
    unlist(blockApply(x, function(bk, na.rm, center) {
        .Call(c_colVars, bk, na.rm, center)
    }, grid=colAutoGrid(x), as.sparse=NA, na.rm=na.rm, center=center, ...))
}

x_rowVars <- function(x, rows=NULL, cols=NULL, na.rm=FALSE, center=NULL, ...,
    useNames=NA)
{
    x_check(x, "Calling SCArray:::x_rowVars() with %s ...")
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    stopifnot(is.null(center) || is.numeric(center))
    k <- x_type(x)
    if (k < 3L)
    {
        x <- x_subset(x, rows, cols)
        if (k == 1L)
            v <- .x_row_vars(x, na.rm, center, ...)
        else
            v <- .x_col_vars(t(x), na.rm, center, ...)
        if (isTRUE(useNames)) names(v) <- rownames(x)
        v
    } else {
        x_msg("Calling DelayedMatrixStats::rowVars() ...")
        callNextMethod()
    }
}

x_colVars <- function(x, rows=NULL, cols=NULL, na.rm=FALSE, center=NULL, ...,
    useNames=NA)
{
    x_check(x, "Calling SCArray:::x_colVars() with %s ...")
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    stopifnot(is.null(center) || is.numeric(center))
    k <- x_type(x)
    if (k < 3L)
    {
        x <- x_subset(x, rows, cols)
        if (k == 1L)
            v <- .x_col_vars(x, na.rm, center, ...)
        else
            v <- .x_row_vars(t(x), na.rm, center, ...)
        if (isTRUE(useNames)) names(v) <- colnames(x)
        v
    } else {
        x_msg("Calling DelayedMatrixStats::colVars() ...")
        callNextMethod()
    }
}

setMethod("rowVars", SMatrix, x_rowVars)
setMethod("colVars", SMatrix, x_colVars)


################

x_rowSds <- function(x, rows=NULL, cols=NULL, na.rm=FALSE, center=NULL, ...,
    useNames=NA)
{
    x_check(x, "Calling SCArray:::x_rowSds() with %s ...")
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    stopifnot(is.null(center) || is.numeric(center))
    k <- x_type(x)
    if (k < 3L)
    {
        x <- x_subset(x, rows, cols)
        if (k == 1L)
            v <- .x_row_vars(x, na.rm, center, ...)
        else
            v <- .x_col_vars(t(x), na.rm, center, ...)
        v <- sqrt(v)
        if (isTRUE(useNames)) names(v) <- rownames(x)
        v
    } else {
        x_msg("Calling DelayedMatrixStats::rowSds() ...")
        callNextMethod()
    }
}

x_colSds <- function(x, rows=NULL, cols=NULL, na.rm=FALSE, center=NULL, ...,
    useNames=NA)
{
    x_check(x, "Calling SCArray:::x_colSds() with %s ...")
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    stopifnot(is.null(center) || is.numeric(center))
    k <- x_type(x)
    if (k < 3L)
    {
        x <- x_subset(x, rows, cols)
        if (k == 1L)
            v <- .x_col_vars(x, na.rm, center, ...)
        else
            v <- .x_row_vars(t(x), na.rm, center, ...)
        v <- sqrt(v)
        if (isTRUE(useNames)) names(v) <- colnames(x)
        v
    } else {
        x_msg("Calling DelayedMatrixStats::colSds() ...")
        callNextMethod()
    }
}

setMethod("rowSds", SMatrix, x_rowSds)
setMethod("colSds", SMatrix, x_colSds)


################

.x_row_w_vars <- function(x, w, na.rm, ...)
{
    # initialize
    stopifnot(is.numeric(w), length(w)==ncol(x))
    if (is.integer(w)) w <- as.double(w)
    .Call(c_init_block)
    # block read
    v <- blockReduce(function(bk, v, w, na.rm) {
        .Call(c_rowWVars, bk, v, w, na.rm)
    }, x, init=matrix(0.0, nrow=nrow(x), ncol=4L), grid=colAutoGrid(x),
        as.sparse=NA, w=w, na.rm=na.rm)
    # finally
    .Call(c_rowWVars_final, v)
}

.x_col_w_vars <- function(x, w, na.rm, ...)
{
    # initialize
    stopifnot(is.numeric(w), length(w)==nrow(x))
    if (is.integer(w)) w <- as.double(w)
    # block read
    unlist(blockApply(x, function(bk, w, na.rm) {
        .Call(c_colWVars, bk, w, na.rm)
    }, grid=colAutoGrid(x), as.sparse=NA, w=w, na.rm=na.rm, ...))
}

x_rowWeightedVars <- function(x, w=NULL, rows=NULL, cols=NULL, na.rm=FALSE,
    ..., useNames=NA)
{
    if (is.null(w))
    {
        x_rowVars(x, rows, cols, na.rm, ..., useNames=useNames)
    } else {
        x_check(x, "Calling SCArray:::x_rowWeightedVars() with %s ...")
        stopifnot(is.logical(na.rm), length(na.rm)==1L)
        k <- x_type(x)
        if (k < 3L)
        {
            x <- x_subset(x, rows, cols)
            if (k == 1L)
                v <- .x_row_w_vars(x, w, na.rm, ...)
            else
                v <- .x_col_w_vars(t(x), w, na.rm, ...)
            if (isTRUE(useNames)) names(v) <- rownames(x)
            v
        } else {
            x_msg("Calling DelayedMatrixStats::rowWeightedVars() ...")
            callNextMethod()
        }
    }
}

x_colWeightedVars <- function(x, w=NULL, rows=NULL, cols=NULL, na.rm=FALSE,
    ..., useNames=NA)
{
    if (is.null(w))
    {
        x_colVars(x, rows, cols, na.rm, ..., useNames=useNames)
    } else {
        x_check(x, "Calling SCArray:::x_colWeightedVars() with %s ...")
        stopifnot(is.logical(na.rm), length(na.rm)==1L)
        k <- x_type(x)
        if (k < 3L)
        {
            x <- x_subset(x, rows, cols)
            if (k == 1L)
                v <- .x_col_w_vars(x, w, na.rm, ...)
            else
                v <- .x_row_w_vars(t(x), w, na.rm, ...)
            if (isTRUE(useNames)) names(v) <- colnames(x)
            v
        } else {
            x_msg("Calling DelayedMatrixStats::colWeightedVars() ...")
            callNextMethod()
        }
    }
}

setMethod("rowWeightedVars", SMatrix, x_rowWeightedVars)
setMethod("colWeightedVars", SMatrix, x_colWeightedVars)


################

x_rowWeightedSds <- function(x, w=NULL, rows=NULL, cols=NULL, na.rm=FALSE,
    ..., useNames=NA)
{
    if (is.null(w))
    {
        x_rowSds(x, rows, cols, na.rm, ..., useNames=useNames)
    } else {
        x_check(x, "Calling SCArray:::x_rowWeightedSds() with %s ...")
        stopifnot(is.logical(na.rm), length(na.rm)==1L)
        k <- x_type(x)
        if (k < 3L)
        {
            x <- x_subset(x, rows, cols)
            if (k == 1L)
                v <- .x_row_w_vars(x, w, na.rm, ...)
            else
                v <- .x_col_w_vars(t(x), w, na.rm, ...)
            v <- sqrt(v)
            if (isTRUE(useNames)) names(v) <- rownames(x)
            v
        } else {
            x_msg("Calling DelayedMatrixStats::rowWeightedSds() ...")
            callNextMethod()
        }
    }
}

x_colWeightedSds <- function(x, w=NULL, rows=NULL, cols=NULL, na.rm=FALSE,
    ..., useNames = NA)
{
    if (is.null(w))
    {
        x_colSds(x, rows, cols, na.rm, ..., useNames=useNames)
    } else {
        x_check(x, "Calling SCArray:::x_colWeightedSds() with %s ...")
        stopifnot(is.numeric(w), length(w)==nrow(x))
        stopifnot(is.logical(na.rm), length(na.rm)==1L)
        k <- x_type(x)
        if (k < 3L)
        {
            x <- x_subset(x, rows, cols)
            if (k == 1L)
                v <- .x_col_w_vars(x, w, na.rm, ...)
            else
                v <- .x_row_w_vars(t(x), w, na.rm, ...)
            v <- sqrt(v)
            if (isTRUE(useNames)) names(v) <- colnames(x)
            v
        } else {
            x_msg("Calling DelayedMatrixStats::colWeightedSds() ...")
            callNextMethod()
        }
    }
}

setMethod("rowWeightedSds", SMatrix, x_rowWeightedSds)
setMethod("colWeightedSds", SMatrix, x_colWeightedSds)


################

.x_row_mins <- function(x, na.rm, ...)
{
    blockReduce(function(bk, v, na.rm) {
        .Call(c_rowMins, bk, v, na.rm)
    }, x, rep(Inf, nrow(x)), grid=colAutoGrid(x), as.sparse=NA, na.rm=na.rm, ...)
}

.x_col_mins <- function(x, na.rm, ...)
{
    unlist(blockApply(x, function(bk, na.rm) {
        .Call(c_colMins, bk, na.rm)
    }, grid=colAutoGrid(x), as.sparse=NA, na.rm=na.rm, ...))
}

x_rowMins <- function(x, rows=NULL, cols=NULL, na.rm=FALSE)
{
    x_check(x, "Calling SCArray:::x_rowMins() with %s ...")
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    switch(x_type(x),
        .x_row_mins(x_subset(x, rows, cols), na.rm),     # 1
        .x_col_mins(t(x_subset(x, rows, cols)), na.rm),  # 2
        {   # 3
            x_msg("Calling DelayedArray::rowMins() ...")
            callNextMethod()
        })
}

x_colMins <- function(x, rows=NULL, cols=NULL, na.rm=FALSE)
{
    x_check(x, "Calling SCArray:::x_colMins() with %s ...")
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    switch(x_type(x),
        .x_col_mins(x_subset(x, rows, cols), na.rm),     # 1
        .x_row_mins(t(x_subset(x, rows, cols)), na.rm),  # 2
        {   # 3
            x_msg("Calling DelayedArray::colMins() ...")
            callNextMethod()
        })
}

setMethod("rowMins", SMatrix, x_rowMins)
setMethod("colMins", SMatrix, x_colMins)


################

.x_row_maxs <- function(x, na.rm, ...)
{
    blockReduce(function(bk, v, na.rm) {
        .Call(c_rowMaxs, bk, v, na.rm)
    }, x, rep(-Inf, nrow(x)), grid=colAutoGrid(x), as.sparse=NA, na.rm=na.rm, ...)
}

.x_col_maxs <- function(x, na.rm, ...)
{
    unlist(blockApply(x, function(bk, na.rm) {
        .Call(c_colMaxs, bk, na.rm)
    }, grid=colAutoGrid(x), as.sparse=NA, na.rm=na.rm, ...))
}

x_rowMaxs <- function(x, rows=NULL, cols=NULL, na.rm=FALSE)
{
    x_check(x, "Calling SCArray:::x_rowMaxs() with %s ...")
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    switch(x_type(x),
        .x_row_maxs(x_subset(x, rows, cols), na.rm),     # 1
        .x_col_maxs(t(x_subset(x, rows, cols)), na.rm),  # 2
        {   # 3
            x_msg("Calling DelayedArray::rowMaxs() ...")
            callNextMethod()
        })
}

x_colMaxs <- function(x, rows=NULL, cols=NULL, na.rm=FALSE)
{
    x_check(x, "Calling SCArray:::x_colMaxs() with %s ...")
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    switch(x_type(x),
        .x_col_maxs(x_subset(x, rows, cols), na.rm),     # 1
        .x_row_maxs(t(x_subset(x, rows, cols)), na.rm),  # 2
        {   # 3
            x_msg("Calling DelayedArray::colMaxs() ...")
            callNextMethod()
        })
}

setMethod("rowMaxs", SMatrix, x_rowMaxs)
setMethod("colMaxs", SMatrix, x_colMaxs)


################

.x_row_ranges <- function(x, na.rm, ...)
{
    blockReduce(function(bk, v, na.rm) {
        .Call(c_rowRanges, bk, v, na.rm)
    }, x, init=t(matrix(c(Inf, -Inf), nrow=2L, ncol=nrow(x))),
        grid=colAutoGrid(x), as.sparse=NA, na.rm=na.rm, ...)
}

.x_col_ranges <- function(x, na.rm, ...)
{
    do.call(rbind, blockApply(x, function(bk, na.rm) {
        .Call(c_colRanges, bk, na.rm)
    }, grid=colAutoGrid(x), as.sparse=NA, na.rm=na.rm, ...))
}

x_rowRanges <- function(x, rows=NULL, cols=NULL, na.rm=FALSE)
{
    x_check(x, "Calling SCArray:::x_rowRanges() with %s ...")
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    switch(x_type(x),
        .x_row_ranges(x_subset(x, rows, cols), na.rm),     # 1
        .x_col_ranges(t(x_subset(x, rows, cols)), na.rm),  # 2
        {   # 3
            x_msg("Calling DelayedArray::rowRanges() ...")
            callNextMethod()
        })
}

x_colRanges <- function(x, rows=NULL, cols=NULL, na.rm=FALSE)
{
    x_check(x, "Calling SCArray:::x_colRanges() with %s ...")
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    switch(x_type(x),
        .x_col_ranges(x_subset(x, rows, cols), na.rm),     # 1
        .x_row_ranges(t(x_subset(x, rows, cols)), na.rm),  # 2
        {   # 3
            x_msg("Calling DelayedArray::colRanges() ...")
            callNextMethod()
        })
}

setMethod("rowRanges", SMatrix, x_rowRanges)
setMethod("colRanges", SMatrix, x_colRanges)


################

.x_row_collapse <- function(x, idx)
{
    if (is.double(idx)) idx <- as.integer(idx)
    ii <- unique(sort(idx))  # NA is removed
    idx <- match(idx, ii)
    x <- x[, ii, drop=FALSE]
    # init
    .Call(c_rowCollapse_init, idx, dim(x))
    on.exit(.Call(c_rowCollapse_done))
    init_val <- rep(NA, nrow(x))
    storage.mode(init_val) <- type(x)
    # block read
    blockReduce(function(bk, v) .Call(c_rowCollapse, bk, v),
        x, init=init_val, grid=colAutoGrid(x))
}

.x_col_collapse <- function(x, idx)
{
    if (is.double(idx)) idx <- as.integer(idx)
    ii <- unique(sort(idx))  # NA is removed
    idx <- match(idx, ii)
    x <- x[ii, , drop=FALSE]
    # init
    .Call(c_colCollapse_init, idx)
    # block read
    blockReduce(function(bk, v) .Call(c_colCollapse, bk, v),
        x, init=vector(type(x), ncol(x)), grid=colAutoGrid(x))
}

x_rowCollapse <- function(x, idxs, rows=NULL, ..., useNames=NA)
{
    x_check(x, "Calling SCArray:::x_rowCollapse() with %s ...")
    stopifnot(is.numeric(idxs), length(idxs) > 0L)
    k <- x_type(x)
    if (k < 3L)
    {
        x <- x_subset(x, rows, NULL)
        if (k == 1L)
            v <- .x_row_collapse(x, idxs)
        else
            v <- .x_col_collapse(t(x), idxs)
        if (isTRUE(useNames)) names(v) <- rownames(x)
        v
    } else {  # 3
        x_msg("Calling DelayedMatrixStats::rowCollapse() ...")
        callNextMethod()
    }
}

x_colCollapse <- function(x, idxs, cols=NULL, ..., useNames=NA)
{
    x_check(x, "Calling SCArray:::x_colCollapse() with %s ...")
    stopifnot(is.numeric(idxs), length(idxs) > 0L)
    k <- x_type(x)
    if (k < 3L)
    {
        x <- x_subset(x, NULL, cols)
        if (k == 1L)
            v <- .x_col_collapse(x, idxs)
        else
            v <- .x_row_collapse(t(x), idxs)
        if (isTRUE(useNames)) names(v) <- colnames(x)
        v
    } else {  # 3
        x_msg("Calling DelayedMatrixStats::colCollapse() ...")
        callNextMethod()
    }
}

setMethod("rowCollapse", SMatrix, x_rowCollapse)
setMethod("colCollapse", SMatrix, x_colCollapse)




