log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

library(data.table)
source("utils/mosaiClassifier/makeSVcalls.R")

probs <- readRDS(snakemake@input[[1]])
probs.biall <- getBiallelicLikelihoods(probs)


fwrite(probs.biall[[1]], file=snakemake@output[["data_table"]], sep="\t", quote=F, row.names=F)
write.table(probs.biall[[2]], file=snakemake@output[["matrix"]], sep="\t", quote=F)

