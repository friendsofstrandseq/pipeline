library(data.table)
library(ggplot2)
library(gridExtra)
source("utils/mosaiClassifier/generateHaploStates.R")
source("utils/test_SVscite/newComputeBreakpointMatrix.R")
source("utils/mosaiClassifier/mosaiClassifier.R")


simulation <- 1 # args[1]
bin.size <- format(100000, scientific = F)#args[2]#
chrom.filt <- "chr1" #args[3]
sce <- fread(paste0("/MMCI/TM/scratch/maryam/SV-SCITE/test-breakpoints-matrix/SV-tree-simulation-2018-01-10/simulation/sce/genome", simulation,"-", bin.size, ".txt"))
sv <- fread(paste0("/MMCI/TM/scratch/maryam/SV-SCITE/test-breakpoints-matrix/SV-tree-simulation-2018-01-10/simulation/phylogeny/variants/genome", simulation, "-", bin.size, ".txt"))
counts <- fread(paste0("zcat /MMCI/TM/scratch/maryam/SV-SCITE/test-breakpoints-matrix/SV-tree-simulation-2018-01-10/simulation/phylogeny/counts/genome", simulation, "-", bin.size, ".txt.gz"))
info <- fread(paste0("/MMCI/TM/scratch/maryam/SV-SCITE/test-breakpoints-matrix/SV-tree-simulation-2018-01-10/simulation/info/genome", simulation, "-", bin.size,".txt"))

hapStatus <- c(hom_ref="1010", hom_del="0000", del_h1="0010", del_h2="1000", hom_dup="2020", dup_h1="2010", dup_h2="1020", hom_inv="0101", inv_h1="0110", inv_h2="1001", inv_dup_h1="1110", inv_dup_h2="1011")
bin.size <- as.numeric(bin.size)

# calling the functions
br.probs <- getBreakpointsLikelihood(sce, counts, chrom.filt, definedHapStatus=T, hapStatus=hapStatus)
sv.br <- getTrueBreakpoints(sv, bin.size, chrom.filt)
sce.br <- getTrueBreakpoints(sce, bin.size, chrom.filt)

br.probs <- addTrueBRtoBRprobs(br.probs, sv.br, sce.br)

# plotting
br.ll.hist <- ggplot(br.probs, aes(x=log(br_ll), fill=breakpoint_type!="no_br"))+geom_histogram(aes(y=..density..))
br.ll.density <- ggplot(br.probs, aes(x=log(br_ll), col=breakpoint_type!="no_br"))+geom_density()
grid.arrange(br.ll.hist, br.ll.density, nrow=2,  ncol=1)



getBreakpointsLikelihood <- function(sce, counts, chrom.filt=NULL, definedHapStatus=F, hapStatus=NULL, haplotypeMode=F){
	# define segs
	num.cells <- length(unique(sce[, cell]))
	segs <- counts[, .(k=.N/num.cells), by=chrom]
	segs <- segs[, cbind(.SD, bps=1:k-1), by=1:nrow(segs)]
	
	# filtering based on chromosome
	if (!is.null(chrom.filt)) {
		counts <- counts[chrom==chrom.filt]
		sce <- sce[chrom==chrom.filt]
	}

	# running mosaiClassifier
	probs <- mosaiClassifierPrepare(counts, info, strand = sce, segs)
	probs <- mosaiClassifierCalcProbs(probs, definedHapStatus=definedHapStatus, hapStatus=hapStatus, haplotypeMode)

	if (!haplotypeMode) {
		# define genotype column
		probs[, genotype:=rep(min(haplotype, sisterHaplotype(haplotype)), .N), by=haplotype]
		# remove the repetitive genotypes
		probs <- probs[haplotype==genotype]
	}

	# sort the probs table based on sample, cell, bin and haplotype
	setkey(probs, sample, cell, chrom, start, end, haplotype)

	# get the number of genotypes
	num.genotypes=length(unique(probs[, genotype]))

	# Add another column including the nb_gt_ll of the next bin (set to -1 for the last bins in chromosomes)
	probs[, next_bin_nb_gt_ll:=data.table::shift(nb_gt_ll, n=num.genotypes, type="lead", fill=-1), 
      by=.(sample, cell, chrom)]

	# compute breakpoint probs
	br.probs <- probs[, .(br_ll=sum(nb_gt_ll*next_bin_nb_gt_ll)), 
                     by=.(sample, cell, chrom, start, end)]

	# remove the start column
	br.probs[, start:=NULL]

	# rename end column to breakpoint
	colnames(br.probs)[which(colnames(br.probs)=="end")]="breakpoint"

	# remove the breakpouint likelihoods for the end bins
	br.probs <- br.probs[br_ll>=0]

	return(br.probs)
}

getTrueBreakpoints <- function(dt, bin.size, chrom.filt=NULL) {
	# filtering based on chromosome
	if (!is.null(chrom.filt)) {
		dt <- dt[chrom==chrom.filt]
	}

	# round start and end coordinates based on the bin.size	
	dt.rounded <- dt	
	dt.rounded[, `:=`(start=round(start/bin.size)*bin.size, end=round(end/bin.size)*bin.size)]

	# add breakpoint column
	dt.rounded <- dt.rounded[, cbind(.SD, breakpoint=c(start, end)),
                           by=.(sample, cell, chrom)]

	# remove start and end columns
	dt.rounded[, `:=`(start=NULL, end=NULL)]

	return(unique(dt.rounded))
}

addTrueBRtoBRprobs <- function(br.probs, sv.br, sce.br) {
	# merging sv and sce breakpoints
	true.br <- merge(sce.br[, -"class"], sv.br, all=T, by=c("sample", "cell", "chrom", "breakpoint"))
	# rename the "SV_type" column to "breakpoint_type"	
	colnames(true.br)[colnames(true.br)=="SV_type"] = "breakpoint_type"
	# set the NA breakpoints to "sce"	
	true.br[is.na(breakpoint_type), breakpoint_type:="sce"]

	# adding true (sce and sv) breakpoints to the br.probs table
	br.probs <- merge(br.probs, true.br, all.x=T, by=c("sample", "cell", "chrom", "breakpoint"))
	# set the NA breakpoints to "no_br"
	br.probs[is.na(breakpoint_type), breakpoint_type:="no_br"]

	return(br.probs)
}
