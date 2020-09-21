# Whoeps, 31st August


library(tidyr)
library(dplyr)
library(data.table)
suppressMessages(library("optparse"))

# Average CN is the CVs weighted by segment length.
weighted_average_CV = function(len_f, avcv_f){
  return(sum((len_f/sum(len_f)) * as.numeric(avcv_f)))
}

# INPUT
#bedfilepath = '/home/hoeps/PhD/projects/huminvs/mosaicatcher/tracks/HG00733.bed'
#querypath = '/home/hoeps/PhD/projects/huminvs/mosaicatcher/tracks/Mapping_Normalization/CN_track_plots/tmp/HG00733_for_CN_survived_bins_sorted.bed'
#invpath = "/home/hoeps/PhD/projects/huminvs/mosaicatcher/tracks/Mapping_Normalization/CN_track_plots/input_invs/HG00733_for_CN.bed"
# PARSE INPUT
option_list = list( 
  make_option(c("-b", "--CNbed"), type="character", default=NULL, 
              help="copynumbers.bed, with 5 columns", metavar="character"),
  make_option(c("-q", "--query"), type="character", default=NULL, 
              help="queries (aka survivors.bed). 5 columns: chr, start, stop, sourceinv, correctly_mapped_reads. No header", metavar="character"),
  make_option(c("-i", "--invfile"), type='character', 
              help='the original inversion file. We want to merge our results with that.'),
  make_option(c("-o", "--outfile"), type="character", default="./output/", 
              help="outputfile", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default=NULL,
	      help="sample name", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
#opt$CNbed ='CN_beds/HG00733_CN.bed'
#opt$query= 'tmp/HG00733_survived_bins_sorted.bed'
#opt$invfile= 'input_invs/HG00733.bed'
#opt$outfile='result/HG00733_CN.txt'
#opt$sample='HG00733'

bedfilepath = opt$CNbed
querypath = opt$query
invpath = opt$invfile
outpath = opt$outfile
sample = opt$sample

# is input file specified?
if (is.null(opt$CNbed)){
  print_help(opt_parser)
  stop("Please specify path to CN bed file", call.=FALSE)
}
if (is.null(opt$sample)){
  print_help(opt_parser)
  stop("Please specify the sample name", call.=FALSE)
}


# TEMP OUTPUT FILE
overlappath = paste0('tmp/',sample,'overlap.bed')

####################################################################################

# Ok back to business.
# assume we have a good bed with 5 columns. Bedmap does the main work for us: get the average CN
# for each entry in the query file, weighted by overlap with that entry in the bigbed. Its a bit
# much to explain now but I promise it works. I tested it. See 
# https://bedops.readthedocs.io/en/latest/content/reference/statistics/bedmap.html for an explanation
#print((paste0('bedmap --echo --wmean <(sort ', querypath,') ', bedfilepath, '> ', overlappath)))
system(paste0('bedmap --echo --wmean ', querypath,' ', bedfilepath, '> ', overlappath))
#Read the survivors.bed as a table so we can later extract the mapping column 
mappings = read.table(querypath, stringsAsFactors = F)
# Read in that file that was just created by bedmap
print("1")
print(overlappath)
overlaps = read.table(overlappath, stringsAsFactors = F)

# Split last column - bedmap formats that a bit weirdly
ff = separate(data = overlaps, col = V5, sep='\\|', into = c("left", "right"))

#rename columns
colnames(ff) = c('chr','start','end','INV', 'mapability','average_CV')


ff$bin_length = (as.numeric(ff$end) - as.numeric(ff$start))
ff = ff %>% group_by(INV) %>% mutate(sum_read_mapability = sum(as.numeric(mapability)),
                                     n_mappable_bins = length(mapability))



# If mapability is not available for some bins (e.g. far outside on the telomere, ignore bins.
ff = ff[ff$average_CV != 'NAN',]

# Take average over CNs 
ff2 = ff %>% group_by(INV) %>% mutate(CN = weighted_average_CV(as.numeric(bin_length), as.numeric(average_CV)))

# cut down to one entry per group, and only the interesting columns
ff_out = ff2[!duplicated(ff2$INV),][,c('INV', 'n_mappable_bins', 'sum_read_mapability', 'CN')]

# Load inversions
invs = read.table(invpath); colnames(invs) = c('chr','start','end','INV')

# Merge CN annotations with INVs. 
annotated_invs = left_join(invs, ff_out, by='INV')

# Now that we have inversion info, we can calculate mapability percent.
annotated_invs$mapability = annotated_invs$sum_read_mapability / (annotated_invs$end - annotated_invs$start)

annotated_invs$sample=sample
annotated_invs_sort = annotated_invs[,c('chr','start','end','sample','INV', 'CN','mapability', 'n_mappable_bins')]

annotated_invs_sort[is.na(annotated_invs_sort$mapability=='NA'),]$n_mappable_bins = 0
annotated_invs_sort[is.na(annotated_invs_sort$mapability=='NA'),]$mapability = 0


# Round CN for cosmetics 
annotated_invs_sort$CN = round(as.numeric(annotated_invs_sort$CN),2)

# Save
write.table(annotated_invs_sort, file=outpath, quote=F, row.names = F, col.names = T)
