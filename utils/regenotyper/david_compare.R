# Compare with David
library(dplyr)
library(stringr)
library(pheatmap)

# Load david's genotypes
david_list = '~/Desktop/desktop_5th_oct/nonred_inversions_n32_genotypes.csv'
dgt = read.table(david_list, header=T, sep=',', stringsAsFactors = F)

# define entries to keep
cols_to_keep = c('seqnames','start','end','width','genoT')

# filter
dgtf = dgt %>%   select(matches(paste(cols_to_keep, collapse="|")))

# rename columns
colnames(dgtf) = str_replace(colnames(dgtf),'genoT_','')


# can we melt together everything?
mygts = melt(callmatrix); 
# make chr name the same (seqnames)
colnames(mygts) = str_replace(colnames(mygts),'chr','seqnames')

# we have to do some more renamings.
colnames(dgtf) = str_replace(colnames(dgtf),'A$','')
colnames(dgtf) = str_replace(colnames(dgtf),'B$','')

davidgts = melt(dgtf, id.vars = c('seqnames','start','end','width'))
colnames(davidgts) = str_replace(colnames(davidgts),'variable','sample')



joined = inner_join(mygts, davidgts, by=c('seqnames','start','end','sample'))
colnames(joined) = c('seqnames','start','end','width','verdict','genotyper','sample','ignore','david')

# rename stuff again
joined$gtyper_simpler = 'complex'
joined[joined$genotyper %in% c('1|1', '1|1_lowconf'),]$gtyper_simpler = 'HOM'
joined[joined$genotyper %in% c('1|0', '1|0_lowconf','0|1', '0|1_lowconf'),]$gtyper_simpler = 'HET'
joined[joined$genotyper %in% c('0|0', '0|0_lowconf'),]$gtyper_simpler = 'REF'
joined[joined$genotyper %in% c('noreads'),]$gtyper_simpler = 'zeroreads'

aa = table(joined$gtyper_simpler, joined$david)
pheatmap(aa)

# Here is a first useful plot
g <- ggplot(joined, aes(david)) + scale_fill_brewer(palette = "Spectral")
g + geom_histogram(aes(fill=gtyper_simpler), 
                   col="black", 
                   size=.1,
                   stat='count') +   # change number of bins
  labs(title="Genotype predictions for 255 inversions in 31 samples") 

# split. 
g <- ggplot(joined[joined$verdict=='pass',], aes(david)) + scale_fill_brewer(palette = "Spectral")
g + geom_histogram(aes(fill=gtyper_simpler), 
                   col="black", 
                   size=.1,
                   stat='count') +   # change number of bins
  labs(title="Passing inversions in 28 matching samples") 

# misorients, to be continued...
miso = cm2[(cm2$ninv/cm$n) > 0.9,]
misojoin = joined[(joined$seqnames %in% miso$chr) & (joined$start %in% miso$start) & (joined$end %in% miso$end),]
g <- ggplot(misojoin, aes(david)) + scale_fill_brewer(palette = "Spectral")
g + geom_histogram(aes(fill=gtyper_simpler), 
                   col="black", 
                   size=.1,
                   stat='count') +   # change number of bins
  labs(title="Passing inversions in 28 matching samples") 


# waehlerwanderung
joined$match = (joined$gtyper_simpler == joined$david)
jmat = cast(joined, seqnames+start+end+width+verdict~sample, value='match')
jmat$match = rowSums(jmat==T)
jmat$nomatch = rowSums(jmat==F)
jmat$matchpct = jmat$match/(jmat$match + jmat$nomatch)
