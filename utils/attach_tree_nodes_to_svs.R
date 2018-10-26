attach_random_tree_nodes_to_svs <- function(genome, phylo.tree, subclonality, seed) {

	# kick out false-dels
	genome <- genome[SV_type!="false_del"]

	# get the root
	root <- setdiff(phylo.tree$edge[,1], phylo.tree$edge[,2])

	# get the set of internal nodes
	internal.nodes <- setdiff(union(phylo.tree$edge[,1], phylo.tree$edge[,2]), union(root, 1:num.cells))

	############ attach SVs to tree nodes
	# get the fraction of the subclonal SVs
	num.svs <- nrow(genome)
	num.subclonal.svs <- floor(subclonality*num.svs)

	# sample a set of subclonal SVs uniformly
	set.seed(seed)
	sampled.svs <- sample(1:nrow(genome), size = num.subclonal.svs)

	# sample a set of internal nodes for attaching SVs
	set.seed(seed)
	sampled.nodes <- sample(internal.nodes, size = num.subclonal.svs)

	# add an extra column indicating clonality to the genome data table
	genome[, is_clonal:=T]
	genome[sampled.svs, is_clonal:=F]

	# add an extra column for the attached tree node id to genome data table (initialized by root for all SVs)
	genome[, sv_tree_node:=root]

	# reassign nodes to the sub_clonal SVs
	genome[!(is_clonal), sv_tree_node:=sort(sampled.nodes)]

	return(list(genome=genome, svs=sampled.svs, nodes=sampled.nodes))
}

create_sv_file <- function(genome, phylo.tree, seed, window.size) {
	sv <- genome
	# get cells
	cells <- phylo.tree$tip.label

	# add sample column
	sv[, sample:=paste0("simulation", seed, "-", format(window.size, scientific = F))]

	# add cell column
	sv <- sv[, cbind(.SD, sapply(.SD[, sv_tree_node], function(x) extract.clade(phy=phylo.tree, node=x)$tip.label))
		   , by = 1:nrow(sv)]
	# naming the new column as cell
	colnames(sv)[length(colnames(sv))] <- "cell"

	# removing extra columns
	sv[, `:=`(nrow=NULL, is_clonal=NULL, sv_tree_node=NULL)]

	return(sv)
}


