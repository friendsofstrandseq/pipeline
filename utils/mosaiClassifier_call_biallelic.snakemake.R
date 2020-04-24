library(data.table)
source("utils/mosaiClassifier/makeSVcalls.R")

probs = readRDS(snakemake@input[["probs"]])
llr   = as.numeric(snakemake@wildcards[["llr"]])

if (snakemake@params[["use_priors"]]) {
probs <- mosaiClassifierPostProcessing(probs)
} else{
probs <- mosaiClassifierPostProcessing(probs, priors=NULL)
}
probs <- forceBiallelic(probs)
tab <- makeSVCallSimple(probs, llr_thr = llr, manual.segs=snakemake@params[["manual_segs"]])

write.table(tab, file = snakemake@output[[1]], sep = "\t", quote=F, row.names = F, col.names = T)
