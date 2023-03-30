#######################################################################
#
# Package name: SCArray
#
# Description:
#     Large-scale single-cell RNA-seq data manipulation with GDS files
#
# Copyright (C) 2023    Xiuwen Zheng (@AbbVie-ComputationalGenomics)
# License: GPL-3
#


#######################################################################
# Generic functions -- scColMeanVar(), scRowMeanVar()

x_num_row_mean_var <- function(x, na.rm=FALSE, useNames=FALSE, ...)
{
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    v <- cbind(rowMeans(x, na.rm=na.rm), rowVars(x, na.rm=na.rm))
    colnames(v) <- c("mean", "var")
    if (isTRUE(useNames)) rownames(v) <- rownames(x)
    v
}

x_num_col_mean_var <- function(x, na.rm=FALSE, useNames=FALSE, ...)
{
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    v <- cbind(colMeans(x, na.rm=na.rm), colVars(x, na.rm=na.rm))
    colnames(v) <- c("mean", "var")
    if (isTRUE(useNames)) rownames(v) <- colnames(x)
    v
}

.x_row_mean_var <- function(x, na.rm, BPPARAM=getAutoBPPARAM())
{
    # block read
    lst <- .parallel_col_reduce(x, BPPARAM,
        Fun = function(bk, v, na.rm) .Call(c_rowVars, bk, v, na.rm, NULL),
        InitFun = .double_nrow3,
        ReduceFun=`+`, na.rm=na.rm)
    # finally
    .Call(c_rowMeanVar_final, lst)
}

.x_col_mean_var <- function(x, na.rm, BPPARAM=getAutoBPPARAM())
{
    # block read
    lst <- .parallel_col_apply(x, BPPARAM,
        Fun = function(bk, na.rm) .Call(c_colMeanVar, bk, na.rm),
        .flatten=FALSE, na.rm=na.rm)
    # finally
    do.call(rbind, lst)
}

x_sc_row_mean_var <- function(x, na.rm=FALSE, useNames=FALSE, ...)
{
    stopifnot(is(x, "SC_GDSMatrix"))
    x_check(x, "Calling SCArray::scRowMeanVar() with %s ...")
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    k <- x_type(x)
    if (k == 1L)
    {
        v <- .x_row_mean_var(x, na.rm, ...)
    } else if (k == 2L)
    {
        v <- .x_col_mean_var(t(x), na.rm, ...)
    } else {
        v <- cbind(rowMeans(x, na.rm=na.rm), rowVars(x, na.rm=na.rm))
    }
    colnames(v) <- c("mean", "var")
    if (isTRUE(useNames)) rownames(v) <- rownames(x)
    v
}

x_sc_col_mean_var <- function(x, na.rm=FALSE, useNames=FALSE, ...)
{
    stopifnot(is(x, "SC_GDSMatrix"))
    x_check(x, "Calling SCArray::scColMeanVar() with %s ...")
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    k <- x_type(x)
    if (k == 1L)
    {
        v <- .x_col_mean_var(x, na.rm, ...)
    } else if (k == 2L)
    {
        v <- .x_row_mean_var(t(x), na.rm, ...)
    } else {
        v <- cbind(colMeans(x, na.rm=na.rm), colVars(x, na.rm=na.rm))
    }
    colnames(v) <- c("mean", "var")
    if (isTRUE(useNames)) rownames(v) <- colnames(x)
    v
}


setGeneric("scRowMeanVar", function(x, na.rm=FALSE, useNames=FALSE, ...)
    standardGeneric("scRowMeanVar"))
setGeneric("scColMeanVar", function(x, na.rm=FALSE, useNames=FALSE, ...)
    standardGeneric("scColMeanVar"))

setMethod("scRowMeanVar", "matrix", x_num_row_mean_var)
setMethod("scColMeanVar", "matrix", x_num_col_mean_var)

setMethod("scRowMeanVar", "Matrix", x_num_row_mean_var)
setMethod("scColMeanVar", "Matrix", x_num_col_mean_var)

setMethod("scRowMeanVar", SMatrix, x_sc_row_mean_var)
setMethod("scColMeanVar", SMatrix, x_sc_col_mean_var)


#######################################################################

x_num_row_nnzero <- function(x, na.counted=NA, ...)
{
    stopifnot(is.logical(na.counted), length(na.counted)==1L)
    x <- x != 0L
    if (isTRUE(na.counted))
    {
        rowSums(x, na.rm=TRUE) + rowSums(is.na(x))
    } else if (isFALSE(na.counted))
    {
        rowSums(x, na.rm=TRUE)
    } else {
        rowSums(x)
    }
}

x_num_col_nnzero <- function(x, na.counted=NA, ...)
{
    stopifnot(is.logical(na.counted), length(na.counted)==1L)
    x <- x != 0L
    if (isTRUE(na.counted))
    {
        colSums(x, na.rm=TRUE) + colSums(is.na(x))
    } else if (isFALSE(na.counted))
    {
        colSums(x, na.rm=TRUE)
    } else {
        colSums(x)
    }
}


.x_row_nnzero <- function(x, na.counted=NA, BPPARAM=getAutoBPPARAM())
{
    .parallel_col_reduce(x, BPPARAM,
        Fun = function(bk, v, na) .Call(c_row_nnzero, bk, v, na),
        InitFun = function(x, ...) integer(nrow(x)),
        ReduceFun=`+`, na=na.counted)
}

.x_col_nnzero <- function(x, na.counted=NA, BPPARAM=getAutoBPPARAM())
{
    .parallel_col_apply(x, BPPARAM,
        Fun = function(bk, na) .Call(c_col_nnzero, bk, na),
        na=na.counted)
}

x_sc_row_nnzero <- function(x, na.counted=NA, ...)
{
    stopifnot(is(x, "SC_GDSMatrix"))
    x_check(x, "Calling SCArray::row_nnzero() with %s ...")
    stopifnot(is.logical(na.counted), length(na.counted)==1L)
    if (x_type(x) == 2L)
    {
        .x_col_nnzero(t(x), na.counted, ...)
    } else {
        .x_row_nnzero(x, na.counted, ...)
    }
}

x_sc_col_nnzero <- function(x, na.counted=NA, ...)
{
    stopifnot(is(x, "SC_GDSMatrix"))
    x_check(x, "Calling SCArray::col_nnzero() with %s ...")
    stopifnot(is.logical(na.counted), length(na.counted)==1L)
    if (x_type(x) == 2L)
    {
        .x_row_nnzero(t(x), na.counted, ...)
    } else {
        .x_col_nnzero(x, na.counted, ...)
    }
}


setGeneric("row_nnzero",
    function(x, na.counted=NA, ...) standardGeneric("row_nnzero"))
setGeneric("col_nnzero",
    function(x, na.counted=NA, ...) standardGeneric("col_nnzero"))

setMethod("row_nnzero", "matrix", x_num_row_nnzero)
setMethod("col_nnzero", "matrix", x_num_col_nnzero)

setMethod("row_nnzero", "Matrix", x_num_row_nnzero)
setMethod("col_nnzero", "Matrix", x_num_col_nnzero)

setMethod("row_nnzero", "DelayedMatrix", x_num_row_nnzero)
setMethod("col_nnzero", "DelayedMatrix", x_num_col_nnzero)

setMethod("row_nnzero", SMatrix, x_sc_row_nnzero)
setMethod("col_nnzero", SMatrix, x_sc_col_nnzero)


#######################################################################

x_runsvd <- function(x, rank, scale=FALSE, approx=TRUE)
{
    # to use crossprod
    if (x_verbose())
        .cat("Using cross product for PCA")
    # fold = 1 for using cross product
    if (approx)
    {
        rv <- BiocSingular::runIrlbaSVD(x, k=rank, nu=rank, nv=rank,
            center=TRUE, scale=scale, deferred=FALSE, fold=1)
    } else {
        rv <- BiocSingular::runExactSVD(x, k=rank, nu=rank, nv=rank,
            center=TRUE, scale=scale, deferred=FALSE, fold=1)
    }
    # output
    out <- list(sdev = rv$d/sqrt(nrow(x) - 1))
    out$rotation <- rv$v
    colnames(out$rotation) <- sprintf("PC%i", seq_len(ncol(out$rotation)))
    out$x <- sweep(rv$u, 2, rv$d, "*")
    colnames(out$x) <- sprintf("PC%i", seq_len(ncol(out$x)))
    out
}


scRunPCA <- function(x, ncomponents=50, ntop=500, subset_row=NULL, scale=FALSE,
    altexp=NULL, name="PCA", exprs_values="logcounts", dimred=NULL,
    n_dimred=NULL, BSPARAM=NULL, BPPARAM=SerialParam(), verbose=TRUE)
{
    # check
    stopifnot(is(x, "SingleCellExperiment") || is(x, "SummarizedExperiment"))
    stopifnot(is.numeric(ncomponents), length(ncomponents)==1L)

    # get the working matrix 'mat'
    y <- x
    if (!is.null(altexp)) y <- altExp(x, altexp)
    transposed <- !is.null(dimred)
    if (!is.null(dimred))
    {
        mat <- reducedDim(y, dimred)
        if (!is.null(n_dimred))
        {
            if (length(n_dimred) == 1L) n_dimred <- seq_len(n_dimred)
            mat <- mat[, n_dimred, drop = FALSE]
        }
    } else {
        mat <- assay(y, exprs_values)
    }

    # save BPPARAM
    oldbp <- getAutoBPPARAM()
    setAutoBPPARAM(BPPARAM)
    on.exit(setAutoBPPARAM(oldbp))
    if (!bpisup(BPPARAM) && !is(BPPARAM, "MulticoreParam"))
    {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM), add=TRUE)
    }
    if (!transposed)
    {
        if (verbose)
            .cat("Select top ", ntop, " features with the highest variances ...")
        cv <- rowVars(mat)
        if (is.null(subset_row))
            subset_row <- sort(head(order(cv, decreasing=TRUE), ntop))
        mat <- mat[subset_row, , drop=FALSE]
        cv <- cv[subset_row]
        if (scale)
        {
            flag <- cv >= 1e-08
            mat <- mat[flag, , drop=FALSE]/sqrt(cv[flag])
            cv <- rep(1, nrow(mat))
        }
        mat <- t(mat)
    } else {
        cv <- colVars(mat)
    }

    # do PCA
    if (verbose)
        cat("Start PCA on the covariance matrix ...\n")
    if (is.null(BSPARAM))
    {
        pca <- x_runsvd(mat, ncomponents, scale=scale)
    } else {
        # use BiocSingular::runPCA instead
        pca <- runPCA(mat, rank=ncomponents, BSPARAM=BSPARAM, BPPARAM=BPPARAM)
    }

    # PCA results
    varExplained <- pca$sdev^2
    percentVar <- varExplained/sum(cv) * 100
    pcs <- pca$x
    rownames(pcs) <- rownames(mat)
    attr(pcs, "varExplained") <- varExplained
    attr(pcs, "percentVar") <- percentVar
    rownames(pca$rotation) <- colnames(mat)
    attr(pcs, "rotation") <- pca$rotation

    # output
    reducedDim(x, name) <- pcs
    x
}

