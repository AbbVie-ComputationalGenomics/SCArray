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

.crossprod_x_row_seq <- function(x)
{
    gd <- rowAutoGrid(x)
    pb <- x_progress(gd)
    if (!is.null(pb)) on.exit(close(pb))
    # block processing
    blockReduce(function(bk, v, pb)
    {
        if (is(bk, "SparseArraySeed")) bk <- as(bk, "CsparseMatrix")
        if (!is.null(pb))
            setTxtProgressBar(pb, currentBlockId())
        .Call(c_add, v, as.matrix(crossprod(bk)))
    }, x, matrix(0, ncol(x), ncol(x)), grid=gd, as.sparse=NA, pb=pb)
}


.crossprod_x_row_parallel <- function(x, sp, bp)
{
    # run
    lst <- bplapply(sp, function(s, x)
    {
        dm <- dim(x)
        if (s[1L]!=1L || s[2L]!=dm[1L])
        {
            # sub rows
            x <- x[seq.int(s[1L], s[2L]), , drop=FALSE]
        }
        gd <- rowAutoGrid(x)
        if (s[3L]==1L && s[4L]==dm[2L])
        {
            # crossprod the whole matrix
            blockReduce(function(bk, v, ii)
            {
                if (is(bk, "SparseArraySeed"))
                    bk <- as(bk, "CsparseMatrix")
                .Call(c_add, v, as.matrix(crossprod(bk)))
            }, x, matrix(0, ncol(x), ncol(x)), grid=gd, as.sparse=NA)
        } else {
            # crossprod the sub matrix
            ii <- seq.int(s[3L], s[4L])  # sub columns
            blockReduce(function(bk, v, ii)
            {
                if (is(bk, "SparseArraySeed"))
                    bk <- as(bk, "CsparseMatrix")
                .Call(c_add, v, as.matrix(crossprod(bk, bk[,ii])))
            }, x, matrix(0, ncol(x), length(ii)), grid=gd, as.sparse=NA, ii=ii)
        }
    }, BPPARAM=bp, x=x)
    # return(list(lst=lst, sp=sp))
    # reduce
    irow <- vapply(sp, `[`, i=3L, 0L)
    iflag <- rep(TRUE, length(irow))
    for (i in unique(irow))
    {
        ii <- which(irow == i)
        k <- ii[1L]
        for (j in ii[-1L])
        {
            lst[[k]] <- .Call(c_add, lst[[k]], lst[[j]])
            iflag[j] <- FALSE
        }
    }
    lst <- lst[iflag]
    # output
    base::Reduce(cbind, lst)
}


# t(x) %*% y, y = x
x_crossprod_x <- function(x, y)
{
    x_check(x, "Calling SCArray:::x_crossprod_x() with %s ...")
    k <- x_type(x)
    if (k == 2L)
    {
        # efficient direction of 'x' is Y (row), result in [ ncol(x), ncol(x) ]
        bp <- getAutoBPPARAM()
        nworker <- .get_num_worker(bp)
        # split rows & columns (not larger than 2GB)
        ngrp <- floor(min(getAutoBlockSize(), 2^31) / (ncol(x)^2*8))
        sp_r <- scNumSplit(nrow(x), min(ngrp, nworker))
        sp_c <- scNumSplit(ncol(x), floor(nworker/length(sp_r)))
        # do
        if (length(sp_r)!=1L || length(sp_c)!=1L)
        {
            x_msg(sprintf(
                "Distributed to %d processes with (%d x %d) sub matrices ...",
                nworker, length(sp_r), length(sp_c)))
            # split list
            sp <- list()
            for (i in sp_r)
                for (j in sp_c) sp <- c(sp, list(c(i, j)))
            # distribute
            .crossprod_x_row_parallel(x, sp, bp)
        } else {
            # sequentially
            .crossprod_x_row_seq(x)
        }
    } else {
        crossprod(as(x, DMatrix))
    }
}


################

# x %*% t(y), y = x
x_tcrossprod_x <- function(x, y)
{
    x_check(x, "Calling SCArray:::x_tcrossprod_x() with %s ...")
    k <- x_type(x)
    if (k == 1L)
    {
        # efficient direction of 'x' is X (col), result in [ nrow(x), nrow(x) ]
        blockReduce(function(bk, v)
        {
            if (is(bk, "SparseArraySeed")) bk <- as(bk, "CsparseMatrix")
            .Call(c_add, v, as.matrix(tcrossprod(bk)))
        }, x, matrix(0.0, nrow(x), nrow(x)), grid=colAutoGrid(x), as.sparse=NA)
    } else {
        tcrossprod(as(x, DMatrix))
    }
}


setMethod("crossprod", c("SC_GDSMatrix", "missing"), x_crossprod_x)
setMethod("tcrossprod", c("SC_GDSMatrix", "missing"), x_tcrossprod_x)

