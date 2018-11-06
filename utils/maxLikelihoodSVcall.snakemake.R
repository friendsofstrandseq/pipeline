log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

library(data.table)
source("utils/mosaiClassifier/makeSVcalls.R")

probs <- readRDS(snakemake@input[[1]])
reg.factor <- snakemake@params[["alt_allele_ll_reg_factor"]]

if (snakemake@wildcards[["sv_call_type"]]=="ML_biallelic_SV_call"){
  probs.biall <- getBiallelicLikelihoods(probs, reg.factor)
  probs <- probs.biall[[1]][nb_gt_ll > 0.5]
  colnames(probs)[which(colnames(probs)=="alt_allele")] <- "geno_name"
}
sv.call <- getMaxLikelihoodSV(probs)

# convert geno names to haplo names since the plotting function does not accept genotype name!!!
#TODO: remove later, plotting function should accept genotype SV calls as well
sv.call[, SV_class:=sapply(SV_class, function(x) gsub("_het", "_h1", x))]

fwrite(sv.call, file=snakemake@output[[1]], sep="\t", quote=F, row.names=F)

