log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

library(data.table)
library(data.tree)
library(ape)
source("utils/attach_tree_nodes_to_svs.R")

seed <- snakemake@wildcards[["seed"]]
window.size <- snakemake@wildcards[["window_size"]]
genome <- fread(snakemake@input[["genome"]])
phylo.tree <- read.tree(snakemake@input[["tree"]])

sv <- create_sv_file(genome, phylo.tree, seed, window.size)

# write the new genome data table (with extra columns)
fwrite(sv, file=snakemake@output[["variants"]], sep="\t")


