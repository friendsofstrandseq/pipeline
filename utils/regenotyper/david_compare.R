# Compare with David
library(dplyr)
library(stringr)
library(pheatmap)

# Load david's genotypes
david_list = '~/PhD/projects/huminvs/mosaicatcher/bed_factory/revision/david_n35_323/nonred_inversions_n35_genotypes.csv'

dgt = read.table(david_list, header=T, sep=',', stringsAsFactors = F)

# define entries to keep
cols_to_keep = c('seqnames','start','end','width','genoT')

# filter
dgtf = dgt %>%   select(matches(paste(cols_to_keep, collapse="|")))

# rename columns
colnames(dgtf) = str_replace(colnames(dgtf),'genoT_','')


# can we melt together everything?
mygts = (melt(callmatrix)); 
# make chr name the same (seqnames)
colnames(mygts) = str_replace(colnames(mygts),'chrom','seqnames')

# we have to do some more renamings.
colnames(dgtf) = str_replace(colnames(dgtf),'A$','')
colnames(dgtf) = str_replace(colnames(dgtf),'B$','')

davidgts = melt(dgtf, id.vars = c('seqnames','start','end','width'))
colnames(davidgts) = str_replace(colnames(davidgts),'variable','sample')



joined = inner_join(mygts, davidgts, by=c('seqnames','start','end','sample'))
colnames(joined) = c('seqnames','start','end','ID','width','valid_bins','genotyper','sample','ignore','david')

# rename stuff again
joined$gtyper_simpler = 'complex'
joined[joined$genotyper %in% c('1|1', '1|1_lowconf'),]$gtyper_simpler = 'HOM'
joined[joined$genotyper %in% c('1|0', '1|0_lowconf','0|1', '0|1_lowconf'),]$gtyper_simpler = 'HET'
joined[joined$genotyper %in% c('0|0', '0|0_lowconf'),]$gtyper_simpler = 'REF'
joined[as.numeric(joined$valid_bins) < 5,]$gtyper_simpler = 'Lowconf'
joined[joined$genotyper %in% c('noreads'),]$gtyper_simpler = 'zeroreads'

aa = table(joined$gtyper_simpler, joined$david)
pheatmap(aa)

# Here is a first useful plot
g <- ggplot(joined, aes(david)) + scale_fill_brewer(palette = "Spectral")
g + geom_histogram(aes(fill=gtyper_simpler), 
                   col="black", 
                   size=.1,
                   stat='count') +   # change number of bins
  labs(title="Genotype predictions") 

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
jmat = cast(joined, seqnames+start+end+width~sample, value='match')
# determine pct match
jmat$match = rowSums(jmat==T)
jmat$nomatch = rowSums(jmat==F)
jmat$matchpct = round(jmat$match/(jmat$match + jmat$nomatch), 3)
jmat[] <- lapply(jmat, as.character)
cm2 = left_join(cm, jmat[,c('seqnames', 'start', 'end', 'match', 'nomatch', 'matchpct')])

# sort and plot
# sort
jmat <- jmat[order(jmat$matchpct),]
jmat$n = 1:dim(jmat)[1]
ggplot(jmat) + geom_point(aes(x=n, y=matchpct)) + labs(title='Inversion match David <-> Genotyper', x='#Inversion',y='Percent matching predictions')

# what happened to where david has mostly 'het'?
# het_mainly_david

dgtf$nhet = rowSums(dgtf=='HET')
het_mainly_david = jmat[jmat$end %in% dgtf[dgtf$nhet > 28,]$end,]
het_mainly_gtyper = callmatrix[callmatrix$end %in% dgtf[dgtf$nhet > 28,]$end,]
jmat_withnhet = left_join(jmat, dgtf[,c('seqnames', 'start', 'end', 'nhet')], by=c('seqnames','start','end'))
jmat_withnhet <- jmat_withnhet[order(jmat_withnhet$matchpct),]
jmat_withnhet$n = 1:dim(jmat_withnhet)[1]

# define jet colormap
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

g = ggplot(jmat_withnhet)  + 
  geom_point(aes(x=n, y=matchpct, color=nhet)) + 
  labs(title='Inversion match David <-> Genotyper, HET', x='#Inversion',y='Percent matching predictions') 
g + scale_fill_gradient(low="blue", high="red")

# plot those
gt2 = melt(het_mainly_gtyper)
g <- ggplot(gt2, aes(verdict))# + scale_fill_brewer(palette = "Spectral")
g + geom_histogram(aes(fill=verdict), 
                   col="black", 
                   size=.1,
                   stat='count')# +   # change number of bins


#misos
#our misos
gtyper_misos = callmatrix[rowSums((callmatrix == '1|1') | (callmatrix=='1|1_lowconf')) >= 28,]
david_misos = dgtf[rowSums(dgtf == 'HOM') > 25,]

gtyper_suspic = callmatrix[callmatrix$end %in% david_misos$end,]

