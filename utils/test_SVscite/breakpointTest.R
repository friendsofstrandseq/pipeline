library(data.table)
library(IRanges)
library(assertthat)
library(dplyr)

source("computeBreakpointMatrix.R")
setwd("../..")
source("utils/mosaiClassifier/mosaiClassifier.R")

hapStatus <- c(hom_ref="1010", hom_del="0000", del_h1="0010", del_h2="1000", hom_dup="2020", dup_h1="2010", dup_h2="1020", hom_inv="0101", inv_h1="0110", inv_h2="1001", inv_dup_h1="1110", inv_dup_h2="1011")

#args = commandArgs(trailingOnly=TRUE)

simulation <- 1 # args[1]
bin.size <- format(100000, scientific = F)#args[2]#
chrom.filt <- "chr1" #args[3]
sce <- fread(paste0("/MMCI/TM/scratch/maryam/SV-SCITE/test-breakpoints-matrix/SV-tree-simulation-2018-01-10/simulation/sce/genome", simulation,"-", bin.size, ".txt"))
sv <- fread(paste0("/MMCI/TM/scratch/maryam/SV-SCITE/test-breakpoints-matrix/SV-tree-simulation-2018-01-10/simulation/phylogeny/variants/genome", simulation, "-", bin.size, ".txt"))
counts <- fread(paste0("zcat /MMCI/TM/scratch/maryam/SV-SCITE/test-breakpoints-matrix/SV-tree-simulation-2018-01-10/simulation/phylogeny/counts/genome", simulation, "-", bin.size, ".txt.gz"))
info <- fread(paste0("/MMCI/TM/scratch/maryam/SV-SCITE/test-breakpoints-matrix/SV-tree-simulation-2018-01-10/simulation/info/genome", simulation, "-", bin.size,".txt"))

if (!is.null(chrom.filt)) {
  counts <- counts[chrom==chrom.filt]
  sv <- sv[chrom==chrom.filt]
  sce <- sce[chrom==chrom.filt]
}

hist(sv[, end-start])
br <- getBreakpointsMatrix(sce, sv, counts, chrom.filt=chrom.filt)
# getting unique rows
br.unq <- unique(br)
#TODO: remove the zero row later on

# define segs
num.cells <- ncol(br)
segs <- counts[, .(k=.N/num.cells), by=chrom]
segs <- segs[, cbind(.SD, bps=1:k-1), by=1:nrow(segs)]

# kick out sce breakpoints and set the strand states to the majority state in chrom/cells
sce <- getMajorityStrandState(sce)

# running mosaiClassifier
d <- mosaiClassifierPrepare(counts, info, strand = sce, segs)
#d.chr1 <- d[chrom == "chr1"]
e <- mosaiClassifierCalcProbs(d.chr1, definedHapStatus=T, hapStatus=hapStatus)

e.wide <- dcast(e, cell+haplotype ~ chrom+start+end, value.var = "nb_gt_ll")
e.wide.consecutive.col.prod <- e.wide[, .(cell, haplotype)]
for(j in 4:ncol(e.wide)){
  #print(j)
  a <- colnames(e.wide)[j-1]
  b <- colnames(e.wide)[j]
  a.split <- strsplit(a,"_")[[1]]
  b.split <- strsplit(b, "_")[[1]]
  # check if the tow bins are consecutive
  if (a.split[1]==b.split[1] && a.split[3]==b.split[2]){
    e.wide.consecutive.col.prod[, eval(a):=e.wide[, .SD[,1]*.SD[,2], .SDcols=c(a,b)]]
  }
}

# summing up the likelihood products in each cell
br.likelihood <- e.wide.consecutive.col.prod[, lapply(.SD, sum), by=cell, .SDcols=3:ncol(e.wide.consecutive.col.prod)]
br.likelihood.mat <- as.matrix(br.likelihood[, 2:ncol(br.likelihood)])
rownames(br.likelihood.mat) <- br.likelihood[, cell]
br.likelihood.mat <- t(br.likelihood.mat)
br.unq

# compare br.likelihoods.mat with br.unq ...
