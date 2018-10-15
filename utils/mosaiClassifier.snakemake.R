log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

library(data.table)
library(assertthat)
source("utils/mosaiClassifier/mosaiClassifier.R")

# The following function converts the bed format to the segs format, which is later used in mosaicatcher
convert_bed_to_segs_format <- function(bed.table, bin.size){
	colnames(bed.table) <- c("chrom", "start", "end")
	segs <- bed.table
	segs[, `:=`(s=ceiling(start/bin.size)-1, e=ceiling(end/bin.size)-1), by=chrom]
	segs[, `:=`(start=NULL, end=NULL)]

	segs <- reshape(segs, direction = "long", varying = c("s", "e"), v.names = "bps", timevar = NULL)
	segs[, id:=NULL]
	setkey(segs, chrom, bps)

	# remove repetitive rows
	segs <- unique(segs)

	segs[, k:=.N, by=chrom]
	setcolorder(segs, c("k", "chrom", "bps"))

	return(segs)
}


# Currently read files from the Snakemake pipeline
counts = fread(paste("zcat",snakemake@input[["counts"]]))
info   = fread(snakemake@input[["info"]])
strand = fread(snakemake@input[["states"]])
segs   = fread(snakemake@input[["bp"]])
bin.size = snakemake@params[["window_size"]]
if (is.null(bin.size)){
	bin.size = as.numeric(snakemake@wildcards[["window_size"]])
}


if (snakemake@params[["manual_segs"]]){
  segs = convert_bed_to_segs_format(segs, bin.size)
}


# DEPERECATED: this version of normalization is no longer used
# is there a normalization file given?
if ("norm" %in% names(snakemake@input) && length(snakemake@input[["norm"]])>0) {
  message("[MosaiClassifier] Read normalization from ", snakemake@input[["norm"]])
  normalization = fread(snakemake@input[["norm"]])
  message("[Warning] Normalization file specified, but this option is no longer available")
} else {
  normalization = NULL
}

# haplotypeMode?
if ("CW" %in% strand$class) {
  haplotypeMode = T
} else {
  haplotypeMode = F
}

d = mosaiClassifierPrepare(counts, info, strand, segs, manual.segs=snakemake@params[["manual_segs"]])
e = mosaiClassifierCalcProbs(d, maximumCN = 4, haplotypeMode = haplotypeMode)

saveRDS(e, file = snakemake@output[[1]])
