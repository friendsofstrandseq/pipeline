source("utils/mosaiClassifier/haploAndGenoName.R")
source("utils/mosaiClassifier/getDispParAndSegType.R")
source("utils/mosaiClassifier/mosaiClassifier.R")


simulateCounts <- function(sv, sce, info, alpha, bin.size, seed){
  # check input data
  # -> test whether sces are disjoint intervals
  test_disjoint_intervals <- function(d) {
    gr <- makeGRangesFromDataFrame(d)
    end(gr) <- end(gr) -1
    return(isDisjoint(gr))
  }
  sce[, assert_that(test_disjoint_intervals(.SD)), by = .(sample, cell)] %>% invisible
  
  
  # renaming sce and sv start and end columns
  colnames(sce)[colnames(sce)=="start"] <- "sce_start"
  colnames(sce)[colnames(sce)=="end"] <- "sce_end"
  colnames(sv)[colnames(sv)=="start"] <- "sv_start"
  colnames(sv)[colnames(sv)=="end"] <- "sv_end"
  
  # merge sce and sv tables
  merged.sce.sv <- merge_sce_sv(sce, sv)
  
  get_range <- function(a, b)
  {
    if (a<=b){
      return(a:b)
    } else {
      return(NULL)
    }
  }
  
  # create counts table initialized from the merged sce and sv data tables
  counts <- merged.sce.sv[,
                 cbind(.SD, bin_start = union(start, get_range(ceiling(start/bin.size),floor(end/bin.size))*bin.size)),
                 by = 1:nrow(merged.sce.sv)]
  counts[, bin_end:=data.table::shift(bin_start, type="lead", fill=max(end)), by=.(cell, chrom)]
  # removing the rows with equal start and end (happens in the cases where there is an sv or sce breakpoint right in a bin boundary)
  counts <- counts[bin_start < bin_end]
  # removing nrow and start and end columns
  counts[, `:=`(nrow=NULL, start=NULL, end=NULL)]
  
  # rename bin_start and bin_end columns
  colnames(counts)[colnames(counts)=="bin_start"] <- "start"
  colnames(counts)[colnames(counts)=="bin_end"] <- "end"
  
  # add a column for SV_type_code
  counts[, SV_type_code:=rep(geno_name_to_code(SV_type), .N),
         by = SV_type]
  
  # add two columns for W and C copy numbers\
  counts <- counts[, cbind(.SD, Wcn=rep(getSegType(class, SV_type_code)[1], .N), 
             Ccn=rep(getSegType(class, SV_type_code)[2], .N)),
         by=.(class, SV_type_code)]
  
  # add NB params to the counts data table
  counts <- merge(counts, info[, .(sample, cell, nb_p, nb_r)], all.X=T, by=c("sample", "cell"))
  
  # add scalar column (set to 1 for now) accounting for mappability factor
  counts[, scalar:=1]
  
  # add expected read counts (W+C) column
  counts[, expected:=((1-nb_p)*nb_r/nb_p)*((end-start)/bin.size)]
  
  # add dispersion parameters separately for W and C
  counts <- add_dispPar(counts)
  
  # generate w and c read counts
  set.seed(seed)
  counts[, `:=`(w=rnbinom(1, size = disp_w, prob = nb_p),
                c=rnbinom(1, size = disp_c, prob = nb_p)),
         by = 1:nrow(counts)]
  
  # sum up the counts in each bin, clean the data table and return the counts
  counts <- counts[, .(start=min(start), end=max(end), class=ifelse(length(unique(class))==1, head(class,1), "?"), c=sum(c), w=sum(w)),
                   by=.(sample, cell, chrom, floor(start/bin.size))]

  # remove the floor column
  counts[, floor:=NULL]
  
  return(counts)
}


merge_sce_sv <- function(sce, sv){
  # creating granges objects of sce and sv
  sce.granges <- GRanges(seqnames = sce[, paste0(sample, "_", chrom, "_", cell)], 
                         ranges = IRanges(start = sce[, sce_start],
                                          end = sce[, sce_end]))
  
  sv.granges <- GRanges(seqnames = sv[, paste0(sample, "_", chrom, "_", cell)], 
                         ranges = IRanges(start = sv[, sv_start],
                                          end = sv[, sv_end]))
  
  # find overlaps between sces and svs
  f <- findOverlaps(sce.granges, sv.granges)
  ovp <- data.table(sce_hits = queryHits(f), sv_hits = subjectHits(f))
  
  # create an expanded table including all overlapping sces and svs
  sce.expanded <- sce[queryHits(f)]
  sv.expanded <- sv[subjectHits(f)]
  cbind.sce.sv <- cbind(sce.expanded, sv.expanded[,.(sv_start, sv_end, SV_type)])
  
  # add missing sce regions (with no SVs)
  cbind.sce.sv <- merge(sce, cbind.sce.sv, all=T)
  cbind.sce.sv[is.na(sv_start),
               `:=`(sv_start=sce_start,
                    sv_end=sce_end,
                    SV_type="hom_ref")]
  
  ### creating the expanded merged table including all sce and sv breakpoints
  # defining start coordinates
  all.breakpoints <- cbind.sce.sv[, cbind(.SD[which.max(sce_end)],
                                          start=sort(setdiff(unique(c(sce_start, sce_end, sv_start, sv_end)), max(sce_end)))),
                                  by = .(sample, cell, chrom)]
  
  # defining end coordinates
  all.breakpoints[, end := data.table::shift(start, type = "lead", fill = max(sce_end)),
                  by = .(sample, chrom, cell)]
  
  # removing useless columns
  all.breakpoints[, `:=`(sce_start=NULL,
                         sce_end=NULL,
                         sv_start=NULL,
                         sv_end=NULL,
                         class=NULL,
                         SV_type=NULL)]
  
  # creating granges object from all breakpoints data table 
  all.granges <- GRanges(seqnames = all.breakpoints[, paste0(sample, "_", chrom, "_", cell)], 
                        ranges = IRanges(start = all.breakpoints[, start], 
                                         end = all.breakpoints[, end]))
  
  # finding overlaps between all breakpoints regions and both sce and sv data tables (in order to fill class and SV-type columns)
  f.sce <- findOverlaps(all.granges, sce.granges, minoverlap = 2)
  f.sv <- findOverlaps(all.granges, sv.granges, minoverlap = 2)
  
  # defining class and SV_type (initialized by ref) columns
  all.breakpoints[, `:=`(class=sce[subjectHits(f.sce), class], 
                         SV_type="hom_ref")]
  
  # correcting SV_type in SV regions
  all.breakpoints[queryHits(f.sv), SV_type:=sv[subjectHits(f.sv), SV_type]]
  
  return(all.breakpoints)
}

