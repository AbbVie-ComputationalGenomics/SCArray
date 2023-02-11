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


# Internal
q_rhdf5 <- quote(suppressPackageStartupMessages(
    require("rhdf5", quietly=TRUE)))
q_HDF5Array <- quote(suppressPackageStartupMessages(
    require("HDF5Array", quietly=TRUE)))


#######################################################################
# Internal functions
#

.cat <- function(...) cat(..., "\n", sep="")

.plural <- function(num) if (num > 1L) "s" else ""

.pretty <- function(x) prettyNum(x, big.mark=",", scientific=FALSE)


# return total # of features, and total # of samples
.gettotnum <- function(gdsfile)
{
    nfeat <- objdesp.gdsn(index.gdsn(gdsfile, "feature.id"))$dim
    if (length(nfeat) != 1L)
        stop("Invalid dimension of 'feature.id'.")
    nsamp <- objdesp.gdsn(index.gdsn(gdsfile, "sample.id"))$dim
    if (length(nsamp) != 1L)
        stop("Invalid dimension of 'sample.id'.")
    c(nfeat, nsamp)
}


# add data.frame or DataFrame to a GDS node
.add_dframe <- function(nd, dat, compress, verbose, prefix="")
{
    for (i in seq_len(NCOL(dat)))
    {
        nm <- paste0(prefix, names(dat)[i])
        nm <- gsub("/", "_", nm, fixed=TRUE)
        obj <- dat[[i]]
        if (is.data.frame(obj) || is(obj, "DataFrame"))
        {
            .add_dframe(nd, obj, compress, verbose, paste0(nm, "."))
        } else if (is.matrix(obj) || is(obj, "Matrix"))
        {
            ss <- colnames(obj)
            if (is.null(ss)) ss <- paste0("c", seq_len(ncol(obj)))
            ss <- make.unique(ss)
            for (j in seq_len(ncol(obj)))
            {
                v <- obj[, j]
                s <- paste0(nm, ".", ss[j])
                if (verbose)
                    .cat("    ", s, "\t[", storage.mode(v), "]")
                add.gdsn(nd, s, v, compress=compress, closezip=TRUE)
            }
        } else {
            if (verbose)
                .cat("    ", nm, "\t[", storage.mode(obj), "]")
            add.gdsn(nd, nm, obj, compress=compress, closezip=TRUE)
        }
    }
    invisible()
}



#######################################################################
# Open a SCArray GDS file
#
scOpen <- function(gdsfn, readonly=TRUE, allow.duplicate=TRUE)
{
    # check
    stopifnot(is.character(gdsfn), length(gdsfn)==1L)
    stopifnot(is.logical(readonly), length(readonly)==1L)
    stopifnot(is.logical(allow.duplicate), length(allow.duplicate)==1L)

    # open the file
    ans <- openfn.gds(gdsfn, readonly=readonly, allow.fork=TRUE,
        allow.duplicate=allow.duplicate, use.abspath=FALSE)

    # check the file format
    a <- get.attr.gdsn(ans$root)
    if (!is.null(a$FileFormat))
    {
        if (identical(a$FileFormat, "SNP_ARRAY"))
            stop("It is a SNP GDS file, please use SNPRelate::snpgdsOpen().")
        if (identical(a$FileFormat, "SEQ_ARRAY"))
            stop("It is a SeqArray GDS file, please use SeqArray::seqOpen().")
        if (!identical(a$FileFormat, "SC_ARRAY"))
            stop("'FileFormat' should be 'SC_ARRAY'")
    }

    # output
    new("SCArrayFileClass", ans)
}


#######################################################################
# Close the SCArray GDS file
#
scClose <- function(gdsfile)
{
    # check
    stopifnot(inherits(gdsfile, "SCArrayFileClass"))
    # close the GDS file
    closefn.gds(gdsfile)
}


#######################################################################
# Get an DelayedArray instance
#
scArray <- function(gdsfile, varname)
{
    # check
    if (is.character(gdsfile))
    {
        f <- openfn.gds(gdsfile, readonly=TRUE, allow.fork=TRUE,
            allow.duplicate=TRUE, use.abspath=FALSE)
        gdsfile <- new("SCArrayFileClass", f)
    } else if (identical(class(gdsfile), "gds.class"))
    {
        # should be exactly gds.class
        gdsfile <- new("SCArrayFileClass", gdsfile)
    }
    if (!inherits(gdsfile, "SCArrayFileClass"))
        stop("'gdsfile' should be a file name, gds.class, or SCArrayFileClass.")
    # new SC_GDSMatrix
    seed <- SCArraySeed(gdsfile, varname)
    as(DelayedArray(seed), "SC_GDSArray")
}


#######################################################################
# Get an SingleCellExperiment/SummarizedExperiment instance
#
scExperiment <- function(gdsfile, sce=TRUE, use.names=TRUE, load.row=TRUE,
    load.col=TRUE)
{
    # check
    if (is.character(gdsfile))
        gdsfile <- scOpen(gdsfile, readonly=TRUE, allow.duplicate=TRUE)
    stopifnot(inherits(gdsfile, "SCArrayFileClass"))
    stopifnot(is.logical(sce), length(sce)==1L)
    stopifnot(is.logical(use.names), length(use.names)==1L)
    stopifnot(is.logical(load.row), length(load.row)==1L)
    stopifnot(is.logical(load.col), length(load.col)==1L)

    # dimnames
    dm <- .gettotnum(gdsfile)
    feat_id <- samp_id <- NULL
    if (use.names)
    {
        feat_id <- read.gdsn(index.gdsn(gdsfile, "feature.id"))
        samp_id <- read.gdsn(index.gdsn(gdsfile, "sample.id"))
    }
    # list all assays
    nm <- ls.gdsn(gdsfile, include.dirs=FALSE)
    x <- vapply(nm, FUN=function(s) {
        identical(objdesp.gdsn(index.gdsn(gdsfile, s))$dim, dm)
    }, FUN.VALUE=TRUE)
    lst <- lapply(nm[x], function(s) {
        m <- scArray(gdsfile, s)
        rownames(m) <- feat_id; colnames(m) <- samp_id
        if (length(dim(m))==2L)
            m <- as(m, SMatrix)
        else
            m <- as(m, SClass)
        m
    })
    names(lst) <- nm[x]

    # load rowData
    rowdat <- NULL
    if (isTRUE(load.row))
    {
        nd <- index.gdsn(gdsfile, "feature.data", silent=TRUE)
        if (!is.null(nd))
        {
            nmlst <- ls.gdsn(nd, include.hidden=FALSE)
            v <- lapply(nmlst, function(nm) read.gdsn(index.gdsn(nd, nm)))
            names(v) <- nmlst
            rowdat <- DataFrame(v, row.names=feat_id)
        }
    }

    # load colData
    coldat <- NULL
    if (isTRUE(load.col))
    {
        nd <- index.gdsn(gdsfile, "sample.data", silent=TRUE)
        if (!is.null(nd))
        {
            nmlst <- ls.gdsn(nd, include.hidden=FALSE)
            v <- lapply(nmlst, function(nm) read.gdsn(index.gdsn(nd, nm)))
            names(v) <- nmlst
            coldat <- DataFrame(v, row.names=samp_id)
        }
    }

    # output
    if (isTRUE(sce))
    {
        # return a SingleCellExperiment object
        if (is.null(coldat))
            SingleCellExperiment(assays=lst, rowData=rowdat)
        else
            SingleCellExperiment(assays=lst, rowData=rowdat, colData=coldat)
    } else {
        # return a SummarizedExperiment object
        if (is.null(coldat))
            SummarizedExperiment(assays=lst, rowData=rowdat)
        else
            SummarizedExperiment(assays=lst, rowData=rowdat, colData=coldat)
    }
}


#######################################################################
# Convert an R object to a single-cell GDS file
# the input R object can be matrix, DelayedMatrix, SummarizedExperiment or
#     SingleCellExperiment
#
scConvGDS <- function(obj, outfn, assay.name=NULL, save.sp=TRUE,
    type=c("float32", "float64", "int32"), compress="LZMA_RA", clean=TRUE,
    verbose=TRUE)
{
    # check
    stopifnot(is.matrix(obj) | inherits(obj, "Matrix") |
        is(obj, "DelayedMatrix") | is(obj, "SingleCellExperiment") |
        is(obj, "SummarizedExperiment"))
    stopifnot(is.character(outfn), length(outfn)==1L)
    stopifnot(is.null(assay.name) || is.character(assay.name))
    stopifnot(is.logical(save.sp), length(save.sp)==1L)
    type <- match.arg(type)
    stopifnot(is.character(compress), length(compress)==1L)
    stopifnot(is.logical(clean), length(clean)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    # row and column names
    rowdat <- coldat <- metadat <- NULL
    nr <- nrow(obj); rownm <- rownames(obj)
    nc <- ncol(obj); colnm <- colnames(obj)
    lst <- list()
    if (is.matrix(obj) || inherits(obj, "Matrix") || is(obj, "DelayedMatrix"))
    {
        assaylst <- list(counts=obj)
    } else if (is(obj, "SingleCellExperiment") | is(obj, "SummarizedExperiment"))
    {
        assaylst <- assays(obj)
        rowdat <- rowData(obj)
        coldat <- colData(obj)
        metadat <- metadata(obj)
    } else {
        stop("Not support the input: ", class(obj)[1L])
    }

    # check row and column names
    if (anyDuplicated(rownm))
    {
        rownm <- make.unique(rownm)
        warning("rownames are not unique, and 'make.unique()' is called",
            immediate.=TRUE)
    }
    if (anyDuplicated(colnm))
    {
        colnm <- make.unique(colnm)
        warning("colnames are not unique, and 'make.unique()' is called",
            immediate.=TRUE)
    }

    # need rownames and colnames
    if (is.null(rownm))
    {
        rownm <- paste0("g", seq_len(nr))
        warning('rownames=NULL, use c("g1", "g2", ...) instead.', immediate.=TRUE)
    }
    if (is.null(colnm))
    {
        colnm <- paste0("c", seq_len(nc))
        warning('colnames=NULL, use c("c1", "c2", ...) instead.', immediate.=TRUE)
    }

    # create a gds file
    if (verbose) .cat("Output: ", outfn)
    outf <- createfn.gds(outfn)
    on.exit(closefn.gds(outf))
    put.attr.gdsn(outf$root, "FileFormat", "SC_ARRAY")
    put.attr.gdsn(outf$root, "FileVersion", "v1.0")
    if (verbose) .cat("Compression: ", compress)

    # add feature and sample IDs
    add.gdsn(outf, "feature.id", rownm, compress=compress, closezip=TRUE)
    add.gdsn(outf, "sample.id", colnm, compress=compress, closezip=TRUE)
    if (verbose) .cat("Dimension: ", nr, " x ", nc)

    # assays
    if (is.null(assay.name)) assay.name <- names(assaylst)
    s <- setdiff(assay.name, names(assaylst))
    assay.name <- intersect(assay.name, names(assaylst))
    assay.name <- setdiff(assay.name, c("feature.id", "sample.id",
        "feature.data", "sample.data", "meta.data"))
    if (verbose)
    {
        .cat("Assay List [", paste(assay.name, collapse=","), "]:")
    }
    if (length(s))
    {
        s <- paste(s, collapse=",")
        warning("No [", s, "] available.", call.=FALSE, immediate.=TRUE)
    }
    # check counts & logcounts
    if (all(c("counts", "logcounts") %in% assay.name))
    {
        warning("'logcounts' should be excluded, since it can be ",
            "generated from 'counts'.", call.=FALSE, immediate.=TRUE)
    }

    # for-each assay
    for (nm in assay.name)
    {
        if (verbose) cat("    ", nm, "\t!", sep="")
        st <- type
        if (save.sp)
        {
            st <- switch(type, float32="sp.real32", float64="sp.real64",
                int32="sp.int32", stop("Invalid 'type': ", type))
        }
        nd <- add.gdsn(outf, nm, valdim=c(nr, 0L), compress=compress,
            storage=st)
        mt <- assaylst[[nm]]
        if (verbose)
            pb <- txtProgressBar(min=0L, max=ncol(mt), width=50L)
        blockApply(mt, function(x) {
            append.gdsn(nd, x)
            if (verbose)
                setTxtProgressBar(pb, start(currentViewport())[2L])
            NULL
        }, grid=colAutoGrid(mt), BPPARAM=NULL)
        readmode.gdsn(nd)
        if (verbose)
        {
            setTxtProgressBar(pb, ncol(mt))
            close(pb)
            cat("    |"); print(nd)
        }
    }

    # add rowData to feature.data
    nfd <- addfolder.gdsn(outf, "feature.data")
    if (!is.null(rowdat))
    {
        if (verbose) .cat("rowData:")
        .add_dframe(nfd, rowdat, compress, verbose)
    }

    # add colData to sample.data
    nfd <- addfolder.gdsn(outf, "sample.data")
    if (!is.null(coldat))
    {
        if (verbose) .cat("colData:")
        .add_dframe(nfd, coldat, compress, verbose)
    }

    # add metadata to meta.data
    nfd <- addfolder.gdsn(outf, "meta.data")
    if (length(metadat) > 0L)
    {
        if (verbose) .cat("metadata:")
        for (i in seq_along(metadat))
        {
            nm <- names(metadat)[i]
            v <- metadat[[i]]
            if (is(v, "DataFrame")) v <- as.data.frame(v)
            if (verbose)
                .cat("    ", nm, "\t[", storage.mode(v), "]")
            add.gdsn(nfd, nm, v, compress=compress, closezip=TRUE)
        }
    }

    # close file
    on.exit()
    closefn.gds(outf)
    if (isTRUE(clean)) cleanup.gds(outfn, verbose=FALSE)
    cat("Done.\n")

    # output
    invisible(normalizePath(outfn))
}


#######################################################################
# Convert CellRanger Market Exchange Format (MEX) files to GDS
#
scMEX2GDS <- function(feature_fn, barcode_fn, mtx_fn, outfn,
    feature_colnm=c("id", "gene", "feature_type"),
    type=c("float32", "float64", "int32"), compress="LZMA_RA", clean=TRUE,
    verbose=TRUE)
{
    # check
    stopifnot(is.character(feature_fn), length(feature_fn)==1L)
    stopifnot(is.character(barcode_fn), length(barcode_fn)==1L)
    stopifnot(is.character(mtx_fn), length(mtx_fn)==1L)
    stopifnot(is.character(outfn), length(outfn)==1L)
    stopifnot(is.character(feature_colnm))
    type <- match.arg(type)
    stopifnot(is.character(compress), length(compress)==1L)
    stopifnot(is.logical(clean), length(clean)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    # load features
    if (anyDuplicated(feature_colnm))
        stop("'feature_colnm' should be unique.")
    if (length(feature_colnm) == 0L)
        feature_colnm <- "id"
    if (verbose) .cat("Load ", sQuote(feature_fn))
    ft <- read.delim(feature_fn, header=FALSE, stringsAsFactors=FALSE)
    if (ncol(ft) > length(feature_colnm))
    {
        feature_colnm <- c(feature_colnm,
            paste0("var", seq.int(length(feature_colnm)+1L,
                length.out=ncol(ft)-length(feature_colnm))))
    }
    colnames(ft) <- feature_colnm[seq_len(ncol(ft))]
    if (!("id" %in% colnames(ft)))
        stop("No id found in the feature file.")
    if (anyDuplicated(ft$id))
        stop("The first column of feature file is used as feature.id, and should be unique.")

    # load barcodes
    if (verbose) .cat("Load ", sQuote(barcode_fn))
    bt <- read.delim(barcode_fn, header=FALSE, stringsAsFactors=FALSE)
    if (anyDuplicated(bt$V1))
        stop("The first column of barcode file is used as sample.id, and should be unique.")

    # load count matrix
    if (verbose) .cat("Load ", sQuote(mtx_fn))
    mt <- readMM(mtx_fn)
    if (verbose)
    {
        nz <- nnzero(mt)
        .cat("    ", nrow(mt), "x", ncol(mt), ", # of nonzeros: ",
            nz, sprintf(" (%.4f%%)", 100*nz/prod(dim(mt))))
    }
    if (nrow(mt) != nrow(ft))
        stop("# of rows should be as the same as feature file.")
    if (ncol(mt) != nrow(bt))
        stop("# of columns should be as the same as # of rows in barcode file.")

    # create a gds file
    if (verbose) .cat("Output: ", outfn)
    outf <- createfn.gds(outfn)
    on.exit(closefn.gds(outf))
    put.attr.gdsn(outf$root, "FileFormat", "SC_ARRAY")
    put.attr.gdsn(outf$root, "FileVersion", "v1.0")
    if (verbose) .cat("Compression: ", compress)

    # add feature and sample IDs
    add.gdsn(outf, "feature.id", ft$id, compress=compress, closezip=TRUE)
    add.gdsn(outf, "sample.id", bt$V1, compress=compress, closezip=TRUE)

    # assays
    if (verbose) .cat("Count matrix ...")
    st <- switch(type, float32="sp.real32", float64="sp.real64",
        int32="sp.int32", stop("Invalid 'type': ", type))
    nd <- add.gdsn(outf, "counts", valdim=c(nrow(mt), 0L), compress=compress,
        storage=st)
    blockApply(mt, function(x) { append.gdsn(nd, x); NULL },
        grid=colAutoGrid(mt), BPPARAM=NULL)
    readmode.gdsn(nd)
    if (verbose) { cat("    |"); print(nd) }

    # add rowData to feature.data
    nfd <- addfolder.gdsn(outf, "feature.data")
    for (i in seq.int(2L, length.out=ncol(ft)-1L))
    {
        add.gdsn(nfd, names(ft)[i], ft[[i]], compress=compress,
            closezip=TRUE)
    }

    # add colData to sample.data
    nfd <- addfolder.gdsn(outf, "sample.data")
    for (i in seq_len(ncol(bt)-1L))
    {
        add.gdsn(nfd, paste0("var", i), bt[[i+1L]], compress=compress,
            closezip=TRUE)
    }

    # close file
    on.exit()
    closefn.gds(outf)
    if (isTRUE(clean)) cleanup.gds(outfn, verbose=FALSE)
    cat("Done.\n")

    # output
    invisible(normalizePath(outfn))
}


#######################################################################
# Convert CellRanger HDF5 files to GDS
#
scHDF2GDS <- function(h5_fn, outfn, group=c("matrix", "mm10"),
    feature_path=character(),
    type=c("float32", "float64", "int32"), compress="LZMA_RA", clean=TRUE,
    verbose=TRUE)
{
    # check
    stopifnot(is.character(h5_fn), length(h5_fn)==1L)
    stopifnot(is.character(outfn), length(outfn)==1L)
    stopifnot(is.character(group), length(group)>=1L)
    stopifnot(is.character(feature_path))
    if (anyDuplicated(feature_path))
        stop("'feature_path' should be unique.")
    type <- match.arg(type)
    stopifnot(is.character(compress), length(compress)==1L)
    stopifnot(is.logical(clean), length(clean)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    # check packages
    if (!eval(q_rhdf5))
        stop("The package 'rhdf5' should be installed.")
    if (!eval(q_HDF5Array))
        stop("The package 'HDF5Array' should be installed.")

    # load HDF5 count matrix
    a <- rhdf5::h5ls(h5_fn, all=TRUE)
    a$group[a$group == "/"] <- ""
    a$path <- substring(paste(a$group, a$name, sep="/"), 2L)
    sep <- ifelse(grepl("/$", group), "", "/")
    if (verbose) .cat("Load ", sQuote(h5_fn))
    g <- intersect(group, a$path)
    if (length(g) == 0L)
        stop("No hdf5 group for sparse matrix.")
    group <- g[1L]
    if (verbose) cat("   ", group)
    mt <- HDF5Array::H5SparseMatrix(h5_fn, group)
    if (verbose)
        .cat("    [", class(mt)[1L], ", ", nrow(mt), " x ", ncol(mt), "]")

    # load barcodes
    nm <- paste(group, "barcodes", sep=sep)
    if (verbose) .cat("    ", nm)
    bt <- rhdf5::h5read(h5_fn, nm)
    if (anyDuplicated(bt))
        stop("'barcodes' should be unique.")
    if (length(bt) != ncol(mt))
        stop("Num. of columns in the count matrix should be # of barcodes.")

    # load features
    if (length(feature_path) == 0L)
    {
        feature_path <- c("genes", "gene_names",
            "features/id", "features/name", "features/feature_type",
            "features/genome")
    }
    test_path <- paste(group, feature_path, sep="/")
    used <- intersect(test_path, a$path)
    if (length(used) == 0L)
        stop("'No valid 'feature_path'.")
    ft <- lapply(used, function(nm) {
        .cat("    ", nm)
        v <- rhdf5::h5read(h5_fn, nm)
        v <- c(v)  # force to regular data type
        if (length(v) != nrow(mt))
            stop(sQuote(nm), " should have ", nrow(mt), " values.")
        v
    })
    names(ft) <- basename(used)
    if (anyDuplicated(ft[[1L]]))
        stop(sQuote(used[1L]), " should be unique.")

    # create a gds file
    if (verbose) .cat("Output: ", outfn)
    outf <- createfn.gds(outfn)
    on.exit(closefn.gds(outf))
    put.attr.gdsn(outf$root, "FileFormat", "SC_ARRAY")
    put.attr.gdsn(outf$root, "FileVersion", "v1.0")
    if (verbose) .cat("Compression: ", compress)

    # add feature and sample IDs
    add.gdsn(outf, "feature.id", ft[[1L]], compress=compress, closezip=TRUE)
    add.gdsn(outf, "sample.id", bt, compress=compress, closezip=TRUE)

    # assays
    if (verbose) .cat("Count matrix ...")
    st <- switch(type, float32="sp.real32", float64="sp.real64",
        int32="sp.int32", stop("Invalid 'type': ", type))
    nd <- add.gdsn(outf, "counts", valdim=c(nrow(mt), 0L), compress=compress,
        storage=st)
    blockApply(mt, function(x) { append.gdsn(nd, x); NULL },
        grid=colAutoGrid(mt), BPPARAM=NULL)
    readmode.gdsn(nd)
    if (verbose) { cat("    |"); print(nd) }

    # add rowData to feature.data
    nfd <- addfolder.gdsn(outf, "feature.data")
    for (i in seq.int(2L, length.out=length(ft)-1L))
    {
        add.gdsn(nfd, names(ft)[i], ft[[i]], compress=compress,
            closezip=TRUE)
    }

    # add colData to sample.data
    addfolder.gdsn(outf, "sample.data")

    # close file
    on.exit()
    closefn.gds(outf)
    if (isTRUE(clean)) cleanup.gds(outfn, verbose=FALSE)
    .cat("Done.")

    # output
    invisible(normalizePath(outfn))
}


#######################################################################
# Wrap the DelayedArray object with SC_GDSArray
#
scObj <- function(obj, verbose=FALSE)
{
    # check
    stopifnot(is.logical(verbose), length(verbose)==1L)
    # do
    if (is(obj, DClass))
    {
        if (is(obj, "DelayedMatrix") && !is(obj, SMatrix))
        {
            obj <- as(obj, SMatrix)
            if (verbose)
                cat("==> SC_GDSMatrix\n")
        } else if (is(obj, DClass) && !is(obj, SClass))
        {
            obj <- as(obj, SClass)
            if (verbose)
                cat("==> SC_GDSArray\n")
        }
        obj
    } else if (is(obj, "SummarizedExperiment"))
    {
        lst <- assays(obj)
        nm <- names(lst)
        for (i in seq_along(lst))
        {
            v <- lst[[i]]
            if (is(v, "DelayedMatrix") && !is(v, SMatrix))
            {
                assay(obj, i) <- as(v, SMatrix)
                if (verbose)
                    .cat(nm[i], " ==> SC_GDSMatrix")
            } else if (is(v, DClass) && !is(v, SClass))
            {
                assay(obj, i) <- as(v, SClass)
                if (verbose)
                    .cat(nm[i], " ==> SC_GDSArray")
            }
        }
        obj
    } else
        stop("obj should be a SummarizedExperiment object.")
}


#######################################################################

scSetMax <- function(x, vmax)
{
    # check
    stopifnot(is(x, "SC_GDSArray"))
    stopifnot(is.numeric(vmax), length(vmax)==1L)
    if (is.na(vmax)) return(x)
    # set a function
    f <- function(a) base::pmin(a, vmax)
    env <- new.env(parent=globalenv())  # need an environment only having vmax
    assign("vmax", vmax, envir=env)
    environment(f) <- env
    # output
    .sc_val(DelayedArray:::stash_DelayedUnaryIsoOpStack(x, f))
}

scSetMin <- function(x, vmin)
{
    # check
    stopifnot(is(x, "SC_GDSArray"))
    stopifnot(is.numeric(vmin), length(vmin)==1L)
    if (is.na(vmin)) return(x)
    # set a function
    f <- function(a) base::pmax(a, vmin)
    env <- new.env(parent=globalenv())  # need an environment only having vmin
    assign("vmin", vmin, envir=env)
    environment(f) <- env
    # output
    .sc_val(DelayedArray:::stash_DelayedUnaryIsoOpStack(x, f))
}

scSetBounds <- function(x, vmin=NaN, vmax=NaN)
{
    # check
    stopifnot(is(x, "SC_GDSArray"))
    stopifnot(is.numeric(vmin), length(vmin)==1L)
    stopifnot(is.numeric(vmax), length(vmax)==1L)
    if (!is.na(vmin) && !is.na(vmax))
    {
        # set a function
        f <- function(a) base::pmin(base::pmax(a, vmin), vmax)
        env <- new.env(parent=globalenv())  # need an environment only having vmin, vmax
        assign("vmin", vmin, envir=env)
        assign("vmax", vmax, envir=env)
        environment(f) <- env
        # output
        .sc_val(DelayedArray:::stash_DelayedUnaryIsoOpStack(x, f))
    } else if (!is.na(vmin))
    {
        scSetMin(x, vmin)  # output
    } else if (!is.na(vmax))
    {
        scSetMax(x, vmax)  # output
    } else
        x
}

scReplaceNA <- function(x, v=0L)
{
    # check
    stopifnot(is(x, "SC_GDSArray"))
    stopifnot(is.numeric(v), length(v)==1L)
    if (is.na(v)) return(x)
    # set a function
    f <- function(a) { a[is.na(a)] <- v; a }
    env <- new.env(parent=globalenv())  # need an environment only having v
    assign("v", v, envir=env)
    environment(f) <- env
    # output
    .sc_val(DelayedArray:::stash_DelayedUnaryIsoOpStack(x, f))
}
