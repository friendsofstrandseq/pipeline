suppressMessages(library(GenomicRanges))
library(data.table)

args = commandArgs(trailingOnly=TRUE)

removeOverlap <- function(sce){
  sce.granges <- GRanges(seqnames = sce[, paste0(sample, "_", chrom, "_", cell)], 
                         ranges = IRanges(start = sce[, start],
                                          end = sce[, end]))
  
  f <- findOverlaps(sce.granges, sce.granges, minoverlap = 2)
  f <- data.table(queryHits(f), subjectHits(f))
  f <- f[V1 < V2]
  
  sce[, max_chrom_coordinate:=rep(max(end), .N), by=.(sample, cell, chrom)]
  d <- unique(sce[f[, V1]])
  chr_cell <- d[, paste0(chrom, "_", cell)]
  sce[paste0(chrom, "_", cell) %in% chr_cell, 
      `:=`(start=min(start), 
           end=max_chrom_coordinate, 
           class=head(class,1)), 
      by=.(sample, chrom, cell)]
  sce <- unique(sce)
  sce[, max_chrom_coordinate:=NULL]
  
  return(sce)
}

sce <- fread(args[1])
sce <- removeOverlap(sce)
fwrite(sce, file=args[1], sep="\t")
