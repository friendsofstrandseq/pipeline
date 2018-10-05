log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

library(data.table)
library(data.tree)
library(ape)
source("utils/attach_tree_nodes_to_svs.R")

seed <- snakemake@wildcards[["seed"]]
num.cells <- snakemake@params[["cell_count"]]
subclonality <- snakemake@params[["subclonality"]]
genome <- fread(snakemake@input[["genome"]])
num.svs <- nrow(genome)

# name the genome columns
colnames(genome) <- c("chrom", "start", "end", "SV_type")

# generate a random phylogenetics tree
set.seed(seed)
phylo.tree <- rtree(num.cells, tip.label = paste0("cell_", 1:num.cells-1))

# remove the edge lengths
phylo.tree$edge.length <- NULL

# attach tree nodes to SVs
genome <- attach_random_tree_nodes_to_svs(genome, phylo.tree, subclonality, seed)

###### writing to the output files
# write the tree to the ouput file in the newick format
write.tree(phylo.tree, file=snakemake@output[["tree"]])

# write the edge sets to the output edges file
write.table(phylo.tree$edge, file=snakemake@output[["edges"]])

# write the new genome data table (with extra columns)
fwrite(genome, file=snakemake@output[["genome"]], sep="\t")


