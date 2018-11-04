log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

library(data.table)
source("utils/mosaiClassifier/makeSVcalls.R")

probs <- readRDS(snakemake@input[[1]])
reg.factor <- as.numeric(snakemake@params[["reg_factor"]])
probs.biall <- getBiallelicLikelihoods(probs, reg.factor)

fwrite(probs.biall[[1]], file=snakemake@output[["data_table"]], sep="\t", quote=F, row.names=F)
write.table(probs.biall[[2]], file=snakemake@output[["matrix"]], sep="\t", quote=F, row.names=F, col.names=F)
write.table(row.names(probs.biall[[2]]), file=snakemake@output[["matrix_rownames"]], quote=F, row.names=F, col.names=F)

