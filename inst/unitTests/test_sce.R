suppressPackageStartupMessages(library(RUnit))
suppressPackageStartupMessages(library(SCArray))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(DelayedMatrixStats))


test_sce_matrix <- function()
{
	# a GDS file for SingleCellExperiment
	fn <- system.file("extdata", "example.gds", package="SCArray")

	sce <- scExperiment(fn)
	mt <- assays(sce)$counts

	cM <- unname(colMeans(mt))
	cM2 <- colMeans2(mt)
	checkEquals(cM, cM2, "counts: colmean")

	rM <- unname(rowMeans(mt))
	rM2 <- rowMeans2(mt)
	checkEquals(rM, rM2, "counts: rowmean")
}
