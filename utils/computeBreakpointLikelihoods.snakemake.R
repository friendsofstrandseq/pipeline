log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

library(data.table)
source("utils/computeBreakpointLikelihoods.R")

sv <- fread(snakemake@input[["variants"]])
sce <- fread(snakemake@input[["sce"]])
counts <- fread(paste("zcat", snakemake@input[["counts"]]))
info <- fread(snakemake@input[["info"]])
alpha <- snakemake@params[["alpha"]]
bin.size <- as.numeric(snakemake@wildcards[["window_size"]])
chrom.filt <- snakemake@wildcards[["chrom"]] #"chr1"

hapStatus <- c(hom_ref="1010", hom_del="0000", del_h1="0010", del_h2="1000", hom_dup="2020", dup_h1="2010", dup_h2="1020", hom_inv="0101", inv_h1="0110", inv_h2="1001", inv_dup_h1="1110", inv_dup_h2="1011")

br.probs <- getBreakpointsLikelihood(sce, counts, chrom.filt, definedHapStatus=T, hapStatus=hapStatus)
sv.br <- getTrueBreakpoints(sv, bin.size, chrom.filt)
sce.br <- getTrueBreakpoints(sce, bin.size, chrom.filt)

br.probs <- addTrueBRtoBRprobs(br.probs, sv.br, sce.br)

gz <- gzfile(snakemake@output[["breakpoint_ll"]], "w")
write.csv(br.probs, gz, row.names=F)
close(gz)
