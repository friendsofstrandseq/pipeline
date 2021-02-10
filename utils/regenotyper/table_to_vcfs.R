# Whoeps, 10th Feb 2021.
# This is the main script for step three of the regenotyper Snakemake. 
# Take an all.txt and turn it into a series of vcfs. 
# Filtering and testing will, i think, be done by another file. 

library(ggplot2)
library(reshape)
library(dplyr)
library(tibble)
library(optparse)

source('clean_genotype_defs.R')
source('vcf_defs.R')


#INPUT INSTRUCTIONS
option_list = list(
  make_option(c("-a", "--alltxt"), type="character", default=NULL,
              help="ArbiGent output (traditionally 'all.txt') to be turned into vcf", metavar="character"),
  make_option(c("-m", "--msc"), type="character", default=NULL,
              help="count debug file from mosaicatcher main run", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default="./outputcorr/",
              help="Outputdir", metavar="character")
)

# Parse input
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

alltxt_file = opt$alltxt
msc_file = opt$msc
outdir = opt$outdir

#alltxt_file = "/home/hoeps/PhD/projects/huminvs/mosaicatcher/analysis/results/U32_freezemerge/sv_probabilities/all.txt"
#msc_file = "/home/hoeps/PhD/projects/huminvs/mosaicatcher/analysis/results/U32_freezemerge/msc.debug"
#outdir = "/home/hoeps/Desktop/arbitrash"


### PARAMETERS ###

# Cutoff for 'lowconf'
save=T
cutoff = 3 
# Second criterion: plus 5
# Complex vs simple: we want complex LLHs to be double the ones
# of simple, and at the same time at least higher by magnitude 5.
bias_add_factor = 5


### FUNCTIONS ###
load_tab <- function(alltxt_file_f){
  # Tiny function to load the input file
  tab = read.table(alltxt_file_f, header=T, stringsAsFactors = F)
  tab = tab %>% mutate(ID = paste0(chrom,'-',start+1,'-INV-',(end-start)+1))
  
  return(tab)
}

simplify_countmatrix <- function(countm){
  # Strongly simplify things
  countm[countm == 'noreads'] <- './.'
  countm[countm == '0000'] <- './.'
  countm[countm == '2200'] <- './.'
  countm[countm == '2110'] <- './.'
  countm[countm == '3100'] <- './.'
  countm[countm == '1101'] <- './1'
  countm[countm == '0100'] <- '1|.'
  countm[countm == '2101'] <- './1'
  countm[countm == '1110'] <- './0'
  countm[countm == '0010'] <- './1'
  countm[countm == '1120'] <- './.'
  countm[countm == '1210'] <- './.'
  countm[countm == '0103'] <- '1|1'
  countm[countm == '1201'] <- './1'
  countm[countm == '0001'] <- './1'
  countm[countm == '1102'] <- './.'
  countm[countm == '1030'] <- '0|0'
  countm[countm == '1200'] <- './.'
  countm[countm == '0301'] <- '1|1'
  countm[countm == '2000'] <- './.'
  countm[countm == '4000'] <- './.'
  countm[countm == '3010'] <- '0|0'
  countm[countm == '0120'] <- '1|0'
  countm[countm == '1100'] <- './.'
  countm[countm == '3001'] <- './1'
  countm[countm == '0030'] <- './.'
  countm[countm == '2020'] <- '0|0'
  countm[countm == '0202'] <- '1|1'
  countm[countm == '1300'] <- './.'
  countm[countm == '1000'] <- './.'
  countm[countm == '2100'] <- './.'
  
  return(countm)
}

load_and_prep_CN <- function(msc_file_f){
  CN = read.table(msc_file, header=1, stringsAsFactors = F)
  CNmerge = as.data.frame(lapply(CN[, c("chrom","start","end","valid_bins")], as.character))
  
  CNmerge = as.tbl(CNmerge)
  CNmerge <- CNmerge %>%  mutate(chrom = as.character(chrom),
                                 start = as.numeric(as.character(start)),
                                 end = as.numeric(as.character(end)),
                                 valid_bins = as.numeric(as.character(valid_bins)))
  return(CNmerge)
}


################
### RUN CODE ###
################

# Load input file
tab = load_tab(alltxt_file)
# Load CN file. It will be used to include 'valid bins' information, which is good to have in the output files.
CNmerge = load_and_prep_CN(msc_file)

# TODO: WHAT IS THIS LINE DOING?
bias_factor = tab$confidence_hard_over_second # This is the first criterion: double 


######### GET 'REPORTED' GENOTYPES ##############
# tabp = [tab]le_[p]rocessed. This has added a 'GT' column that is the GT result
# of choice. And this GT of choice depends on: bias(add)factor, cutoff and wether
# or not we want 'lowconf' label included.

# tabp = Complex calls allowed, lowconf label added
tabp = add_gts_revisited_lowconf(tab, bias_factor, bias_add_factor, cutoff)
# tabp2 =  Complex calls allowed, lowconf label nope, LLHs printed
tabp2 = add_long_gts_revisited(tab, bias_factor, bias_add_factor, cutoff)
# tabp3 = Complex calls allowed, lowconf label nope
tabp3 = add_gts_revisited(tab, bias_factor, bias_add_factor, cutoff)


# Merge valid bin information into tabp's
tabp <- full_join(tabp, CNmerge, by = c("chrom","start","end"))
tabp3 <- full_join(tabp3, CNmerge, by = c("chrom","start","end"))

####################################################################################

# Cast table tabp3 into vcf-like matrix
callmatrix = cast(unique(tabp3), chrom+start+end+ID+len+valid_bins~sample, value='GT')

####################################################################################

# Name shortening, a bit of reformatting
cms = callmatrix
samplenames = colnames(cms)[7:length(colnames(cms))]
cms_gts = cms[,samplenames]
cms_gts[] <- lapply(cms_gts, as.character)
# Complex calls to simple ones
countm = simplify_countmatrix(cms_gts)
# Bind description and GTs back together
cms_full = cbind(cms[,1:6], cms_gts)

####################################################################################

# Tabp: this is the 'normal' one. With lowconf label
callmatrix = cast(unique(tabp), chrom+start+end+ID+len+valid_bins~sample, value='GT')

# Tabp2: this is the one with LLH info included
callmatrix_detail = cast(tabp2, chrom+start+end+ID+len~sample, value='GTL')
callmatrix_detail = callmatrix_detail[ , colSums(is.na(callmatrix_detail)) == 0]

# Sidequest: find hom invs. A bit messy but we keep it for now. 
callmatrix_hom_lab = callmatrix
callmatrix_hom_lab$nhom = rowSums(callmatrix_hom_lab=='1|1')
cm_hom = callmatrix_hom_lab[callmatrix_hom_lab$nhom > (dim(callmatrix_hom_lab)[2] - 5)*0.8,]
hom_ids = cm_hom$ID
cm_detail_hom = callmatrix_detail[callmatrix_detail$ID %in% hom_ids,]

################## MAKE VCFS
# The detailed one from tabp2, including LLH info.
vcf = vcfify_callmatrix_detail(callmatrix_detail)
# Like above, but filtered for misos
vcf_miso = vcfify_callmatrix_detail(cm_detail_hom)
# And this one is based on simplified tabp3.
vcf_limix = vcfify_callmatrix_simple_for_limix(cms_full)

################## Make at least one plot TODO

#ggplot(tabp3) + geom_point(aes(x=log1p(confidence_hard_over_second), y=log1p(confidence_nobias_over_hard), col=GT)) 
  

################# SAVE ALL THESE DIFFERENT THINGS.
if (save==T){
  
  # Prep directory
  dir.create(outdir)
  
  # Paths, paths, paths. 
  callmatrix_file = 'res.csv'
  callmatrix_file_detail = 'res_detail.csv'
  vcffile_all = 'res_all.vcf'
  vcffile_miso = 'res_miso.vcf'
  vcffile_limix = 'res_verysimple.vcf'
  
  # Save simple callmatrix
  write.table(callmatrix, file=file.path(outdir,callmatrix_file), quote = F, col.names = T, row.names=F, sep='\t')
  # Same detailed callmatrix
  write.table(callmatrix_detail, file=file.path(outdir, callmatrix_file_detail), quote = F, col.names = T, row.names=F, sep='\t')
  # Save vcf_all
  outvcffile = file.path(outdir, vcffile_all)
  writeLines(vcf[[1]], file(outvcffile))
  write.table(vcf[[2]], file=outvcffile, quote=F, col.names=F, row.names=F, sep='\t', append = T)
  system(paste0('bcftools sort -O v -o ', outvcffile, ' ', outvcffile))
  
  # Save vcf_miso
  outvcffile_miso = file.path(outdir, vcffile_miso)
  writeLines(vcf_miso[[1]], file(outvcffile_miso))
  write.table(vcf_miso[[2]], file=outvcffile_miso, quote=F, col.names=F, row.names=F, sep='\t', append = T)
  system(paste0('bcftools sort -O v -o ', outvcffile_miso, ' ', outvcffile_miso))
  
  # Save vcf_limix
  outvcffile_limix = file.path(outdir, vcffile_limix)
  writeLines(vcf_limix[[1]], file(outvcffile_limix))
  write.table(vcf_limix[[2]], file=outvcffile_limix, quote=F, col.names=F, row.names=F, sep='\t', append = T)
  system(paste0('bcftools sort -O v -o ', outvcffile_limix, ' ', outvcffile_limix))
  
  # Save only complex stuff Complex exploration
  tabcomp = tab[tab$GTs == 'complex',]
  indiv_invs <- unique( tabcomp[ , 1:3 ] )
  library(dplyr)
  aa = left_join(indiv_invs,callmatrix_detail)
  bb = aa[
    with(aa, order(chrom, start)),
    ]
}


