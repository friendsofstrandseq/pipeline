library(data.table)
library(GenomicRanges)
library(ggpubr)

biall.probs <- fread("biallelic-likelihood-table.data") #fread(snakemake@input[["biall_probs_table"]])
biall.probs.mat <- fread("biallelic-likelihood-matrix.data") #fread(snakemake@input[["biall_probs_matrix"]])
biall.probs.mat.rownames <- fread("biallelic-likelihood-matrix-rownames.data", header=F)[, V1] #fread(snakemake@input[["biall_probs_matrix_rownames"]])
true.svs <- fread("../../../simulation/phylogeny/genome/genome2.tsv")#fread(snakemake@input[["genome"]])
num.cells <- ncol(biall.probs.mat)

correctSVnames <- function(true.svs){
### FIXME: make the namings fot simulated anaf true SVs consistent
  true.svs[SV_type=="inv_dup", SV_type:="idup_het"]
  true.svs[SV_type=="hom_inv", SV_type:="inv_hom"]
  true.svs[SV_type=="hom_dup", SV_type:="dup_hom"]
  true.svs[SV_type=="het_dup", SV_type:="dup_het"]
  true.svs[SV_type=="het_inv", SV_type:="inv_het"]
  true.svs[SV_type=="het_del", SV_type:="del_het"]
  true.svs[SV_type=="hom_del", SV_type:="del_hom"]
}

mapDetectedToTrueSVs <- function(biall.probs, true.svs, min.seg.cover=0.8)
{
  # get granges
  true.svs.gr <- GRanges(unique(true.svs[,.(chrom, start, end, SV_type, is_clonal, sv_tree_node)]))
  
  segs.gr <- GRanges(unique(biall.probs[,.(chrom, start, end, alt_allele)]))
  
  # finding overlaps between the true and the detected SVs
  ovp <- findOverlaps(segs.gr, true.svs.gr)
  
  # create a data.table of the overlapping segments
  overlap <- data.table(chrom=as.character(seqnames((segs.gr[queryHits(ovp)]))),
			pred.start=start(segs.gr[queryHits(ovp)]),
                        pred.end=end(segs.gr[queryHits(ovp)]),
			pred.sv=segs.gr$alt_allele[queryHits(ovp)],
                        true.start=start(true.svs.gr[subjectHits(ovp)]), 
                        true.end=end(true.svs.gr[subjectHits(ovp)]),
			true.sv=true.svs.gr$SV_type[subjectHits(ovp)],
			is.clonal=true.svs.gr$is_clonal[subjectHits(ovp)],
			sv.tree.node=true.svs.gr$sv_tree_node[subjectHits(ovp)])

  # add a column for overlap length
  overlap[, overlap.length:=min(pred.end, true.end)-max(pred.start, true.start),
            by=1:nrow(overlap)]

  # keep overlaps that have the same pred and true SV and cover at least min.seg.cover fraction of the segment
  overlap <- overlap[overlap.length >= min.seg.cover*(pred.end-pred.start) & pred.sv==true.sv]

  # create a data table mapping pred svs to true svs\
  d <- data.table(pred.sv=overlap[, paste0(chrom, "_", pred.start, "_", pred.end, "_", pred.sv)],
                  true.sv=overlap[, paste0(chrom, "_", true.start, "_", true.end, "_", true.sv)],
		  is.clonal=overlap[,is.clonal],
		  sv.tree.node=overlap[,sv.tree.node])

  return(d)
}


#biall.probs.mat <- fread("biallelic-likelihood-matrix.data")
correctSVnames(true.svs)
overlap <- mapDetectedToTrueSVs(biall.probs, true.svs)

#biall.probs.mat[, is_true_sv := V1 %in% overlap[, pred.sv]]

biall.probs[, pred.sv := paste0(chrom, "_", start, "_", end, "_", alt_allele)]
biall.probs <- merge(biall.probs, overlap[, .(pred.sv, is.clonal)], all.x=T, by="pred.sv")
biall.probs[, sv_status:=ifelse(is.na(is.clonal),"no_sv", ifelse(is.clonal, "clonal_sv", "sub_clonal_sv"))]
biall.probs <- biall.probs[, -"is.clonal"]

biall.probs[, agg_log_ll:=sum(log(nb_gt_ll)), by=.(sample, chrom, start, end)]

#biall.probs[, is_true_sv := paste0(chrom, "_", start, "_", end, "_", alt_allele) %in% overlap[, pred.sv]]
#setkey(biall.probs.mat, is_true_sv)
true.sv.heatmap <- ggplot(biall.probs[sv_status!="no_sv"], aes(cell,paste0(chrom, "_", start, "_", end, "_", alt_allele))) + geom_tile(aes(fill=nb_gt_ll), color="white") + scale_fill_gradient(low="white", high="steelblue")
false.sv.heatmap <- ggplot(biall.probs[sv_status=="no_sv"], aes(cell,paste0(chrom, "_", start, "_", end, "_", alt_allele))) + geom_tile(aes(fill=nb_gt_ll), color="white") + scale_fill_gradient(low="white", high="steelblue")
ggarrange(true.sv.heatmap, false.sv.heatmap, nrow=1, ncol=2)

min.alt.ll.cutoff <- c(0.5, 0.6, 0.7, 0.8, 0.9)
for (i in min.alt.ll.cutoff) {
  biall.probs[, paste0("num_cells_passed_cutoff_", i):=nrow(.SD[nb_gt_ll > i]), by=.(sample, chrom, start, end)]
}

biall.probs.agg <- biall.probs[, head(.SD,1), by=.(sample, chrom, start, end)]

plt.cutoff0.5 <- ggplot(biall.probs.agg, aes(x=num_cells_passed_cutoff_0.5, y = ..density.., col=sv_status))+geom_freqpoly()+ggtitle("num_cells_passed_cutoff_0.5")
plt.cutoff0.6 <- ggplot(biall.probs.agg, aes(x=num_cells_passed_cutoff_0.6, y = ..density.., col=sv_status))+geom_freqpoly()+ggtitle("num_cells_passed_cutoff_0.6")
plt.cutoff0.7 <- ggplot(biall.probs.agg, aes(x=num_cells_passed_cutoff_0.7, y = ..density.., col=sv_status))+geom_freqpoly()+ggtitle("num_cells_passed_cutoff_0.7")
plt.cutoff0.8 <- ggplot(biall.probs.agg, aes(x=num_cells_passed_cutoff_0.8, y = ..density.., col=sv_status))+geom_freqpoly()+ggtitle("num_cells_passed_cutoff_0.8")
plt.cutoff0.9 <- ggplot(biall.probs.agg, aes(x=num_cells_passed_cutoff_0.9, y = ..density.., col=sv_status))+geom_freqpoly()+ggtitle("num_cells_passed_cutoff_0.9")
plt.agg.log.ll <- ggplot(biall.probs.agg, aes(x=agg_log_ll, y = ..density.., col=sv_status))+geom_freqpoly()+ggtitle("agg_log_ll")
plt.agg.ll <- ggplot(biall.probs.agg, aes(x=exp(agg_log_ll), y = ..density.., col=sv_status))+geom_freqpoly()+ggtitle("agg_ll")

plt.save <- ggarrange(plt.cutoff0.5, plt.cutoff0.6, plt.cutoff0.7, plt.cutoff0.8, plt.cutoff0.9, plt.agg.log.ll, nrow=3, ncol=2)
ggsave(plt.save+theme(plot.title=element_text(hjust=0.5)), file="comparative_plots.pdf")

# keep only the segments where the fraction of cells with alt allele prob at least 0.9 are between 0.2 and 0.8
biall.probs.agg <- biall.probs.agg[num_cells_passed_cutoff_0.9 > 0.2*num.cells & num_cells_passed_cutoff_0.9 < 0.8*num.cells]

# make the matrix small
biall.probs.mat <- biall.probs.mat[biall.probs.mat.rownames %in% biall.probs.agg[, paste0(chrom, "_", start, "_", end, "_", alt_allele)]]
biall.probs.mat.rownames <- biall.probs.mat.rownames[biall.probs.mat.rownames %in% biall.probs.agg[, paste0(chrom, "_", start, "_", end, "_", alt_allele)]]
### FIXME: save it as an RData object also in the previous rule where it was created
### remove the following four lines after this fix
#rnames <- biall.probs.mat[, V1]
#biall.probs.mat <- biall.probs.mat[, -"V1"]
#colnames(biall.probs.mat) <- paste0("cell_",0:(num.cells-1))
#row.names(biall.probs.mat) <- rnames

write.table(biall.probs.mat, file="biallelic-likelihood-selected-matrix.data", sep="\t", quote=F, row.names=F, col.names=F) # write.table(biall.probs.mat, file=snakemake@output[["selected_mat"]], sep="\t", quote=F)
write.table(biall.probs.mat.rownames, file="biallelic-likelihood-selected-matrix-rownames.data", quote=F, row.names=F, col.names=F)
fwrite(overlap, file="predToTrueSVmap.data", sep="\t", quote=F, row.names=F)


#test0.5 <- biall.probs.agg[num_cells_passed_cutoff_0.5>20 & num_cells_passed_cutoff_0.5<80]
#test0.6 <- biall.probs.agg[num_cells_passed_cutoff_0.6>20 & num_cells_passed_cutoff_0.6<80]
#test0.7 <- biall.probs.agg[num_cells_passed_cutoff_0.7>20 & num_cells_passed_cutoff_0.7<80]
#test0.8 <- biall.probs.agg[num_cells_passed_cutoff_0.8>20 & num_cells_passed_cutoff_0.8<80]
#test0.9 <- biall.probs.agg[num_cells_passed_cutoff_0.9>20 & num_cells_passed_cutoff_0.9<80]
#table(test0.5[, sv_status])
#table(test0.6[, sv_status])
#table(test0.7[, sv_status])
#table(test0.8[, sv_status])
#table(test0.9[, sv_status])

