# define a function for computing a matrix breakpoints/cells with values 1 for sv breakpoints, 2 for sce breakpoints, and 0 otherwise
getBreakpointsMatrix <- function(sce, sv, counts, chrom.filt=NULL){
  bin.size <- median(counts[, end-start])
  cell.names <- unique(counts[, cell])
  num.cells <- length(cell.names)
  num.bins <- nrow(counts) / num.cells
  bins <- counts[, head(.SD, 1), by=.(chrom, start, end)][,.(chrom, start, end)]
  bin.names <- paste0(bins[, chrom], "_", bins[, start], "_", bins[, end])
  chrom.bin.counts <- bins[, .N, by=chrom]
  # test whether chroms are sorted
  if (is.null(chrom.filt)) {  
    assert_that(all(chrom.bin.counts[, chrom] == c(paste0("chr", c(1:22, "X", "Y"))))) %>% invisible()
  }
  chrom.bin.cumsum.counts <- cumsum(chrom.bin.counts[, N])
  chrom.bin.counts[, cumsum_counts:=cumsum(N)]
  chrom.bin.counts[, chrom_start:=data.table::shift(cumsum_counts, n=1, fill=0, type="lag")]
  # create a zero matrix
  m <- matrix(0L, nrow = num.bins, ncol = num.cells, dimnames = list(bin.names, cell.names))
  
  addBinIndex <- function(df, chrom.bin.counts){
    df <- merge(df, chrom.bin.counts[,.(chrom, chrom_start, cumsum_counts)], by="chrom")
    df[, `:=`(start_bin_idx = chrom_start+round(start/bin.size), end_bin_idx = chrom_start+round(end/bin.size))]
  }
  
  sv <- addBinIndex(sv, chrom.bin.counts)
  if(!is.null(sce)){
    sce <- addBinIndex(sce, chrom.bin.counts)
  }
  
  # try to vectorize and optimize this function later
  fillBreakPointMatrix <- function(df, matrix, fill){
    for (i in 1:nrow(df)){
      for (b in unique(df[i, c(start_bin_idx, end_bin_idx)])){
        if (!b %in% c(0, df[i, cumsum_counts])){
          matrix[b, df[i, cell]] <- fill
        }
      }
    }
    return(matrix)
  }
  
  # add SV breakpoints
  m <- fillBreakPointMatrix(sv, matrix=m, fill=1)
  # add SCE breakpoints
  m <- fillBreakPointMatrix(sce, matrix=m, fill=2)
  
  # TODO (check out later): there are a lot of SVs happening inside one single bin
  
  return(m)
}

# The following function computes the majority strand state in each chromosome and single-cell
getMajorityStrandState <- function(sce)
{
  sce[, `:=`(majority_class_width=max(end-start), chrom_start=min(start), chrom_end=max(end)), by=.(chrom, cell)]
  sce <- sce[end-start==majority_class_width]
  sce[, `:=`(chrom_start=NULL, chrom_end=NULL, majority_class_width=NULL)]
  return(sce)
}

