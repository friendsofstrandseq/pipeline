# Whoeps, 31 May 2020
# This is a perfectly input-output-oriented script. 

# Input: - probabilities.R
#        - bed file with groups
#        - outdir
#        - CN_mapabilty per segment 
#        
# Output: - table with classifications, respective CN and mapability
#         - dumbbellplot
#         - beeswarmplots
print('Initialising Regenotyper...')
# supposed to suppress warnings, but not working
oldw <- getOption("warn")
options(warn = -1)

#!/usr/bin/env Rscript
suppressMessages(library("optparse"))
suppressMessages(library("tidyr"))
suppressMessages(library("stringr"))
suppressMessages(source("probability_helpers.R"))
suppressMessages(source("regenotype_helpers.R"))
suppressMessages(source("bulk_helpers.R"))
# Slighlty obscure package to handle numbers smaller than 10^-300 (LLHs can get small...)
suppressMessages(library("Brobdingnag"))


# INPUT INSTRUCTIONS
option_list = list( 
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="probabilities.R, produced by mosaicatcher", metavar="character"),
  make_option(c("-b", "--bed"), type="character", default=NULL, 
              help="a bed file specifiying labels/groups for the segments.", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default="./outputcorr/", 
              help="output dir name [default= %default]", metavar="character"),
  make_option(c("-c", "--cn_map"), type="character", default=NULL, 
              help="average copy numbers and mapability for all segments in given bed file", metavar="character")
); 


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);



# Ok lets go! Calm the minds of impatient humans first of all.
print('Processing and summarizing information, making plots')
print('This can take a few minutes.')

### This block is for debug mode. ###
#p_link = '/home/hoeps/PhD/projects/huminvs/mosaicatcher/sv_probabilities/HG00733/100000_fixed_norm.selected_j0.1_s0.5/probabilities.Rdata'
#p_link ='../../results/HG733_norm2/probabilities.Rdata'
#p_link = '../../results/HG00733_wolfy2/probabilities.Rdata'
#outdir_raw = '../../results/HG00733_wolfy2/restrash'
#labels_link = '../../results/HG733_hgsvc/naming_HG733_sd.bed'
#outdir_raw = '../../results/HG00733_diagnose/resbulk2/'
#labels_link = '/home/hoeps/PhD/projects/huminvs/mosaicatcher/analysis/results/HG733_postcorr/naming_HG733.bed'
#####################################
#p_link ='../../results/U24_peter/HG00733_aug31_redo/probabilities.Rdata'
#p_link = '../../results/U24_peter/100000_fixed_norm.selected_j0.1_s0.1/probabilities.Rdata'
#labels_link = '../../results/U24_peter/naming/HG00733_peterash_sort.bed'
#outdir_raw = '../../results/deleteme4/'
#outdir_raw = '../../results/deletemealso3/'
#labels_link = '/home/hoeps/PhD/projects/huminvs/mosaicatcher/analysis/results/HG01114_first/easy_naming.bed'
#p_link = '../../results/U24_peter/sv_probabilities/HWWKWAFXY_HG03683x01_19s004569-1-1/100000_fixed_norm.selected_j0.1_s0.5/probabilities.Rdata'
#labels_link = '/home/hoeps/PhD/projects/huminvs/mosaicatcher/bed_factory/audano_3/wgot/done_withnames/newnames/HG00733.bed' 
#outdir_raw = '../../results/U24_peter/sv_probabilities/HWWKWAFXY_HG03683x01_19s004569-1-1/res/'
#p_link = "/home/hoeps/PhD/projects/huminvs/mosaicatcher/analysis/results/david_hack_wolfy/sv_probabilities/HWWKWAFXY_HG03683x01_19s004569-1-1/100000_fixed_norm.selected_j0.1_s0.5/probabilities.Rdata"
#pp = readRDS(p_link)
#p_link = "/home/hoeps/PhD/projects/huminvs/mosaicatcher/analysis/results/U24_peter/sv_probabilities/HWWKWAFXY_HG03683x01_19s004569-1-1/100000_fixed_norm.selected_j0.1_s0.5/probabilities.Rdata"
#p_link = "/home/hoeps/PhD/projects/huminvs/mosaicatcher/analysis/results/sv_probabilities_U24_David_calls_hackathon/HWWKWAFXY_HG03683x01_19s004569-1-1/100000_fixed_norm.selected_j0.1_s0.1/probabilities.Rdata"
#p_link = "/home/hoeps/PhD/projects/huminvs/mosaicatcher/analysis/results/hack_32/sv_probabilities/HG00733/100000_fixed_norm.selected_j0.1_s0.1/probabilities.Rdata"
#p_link = "/home/hoeps/PhD/projects/huminvs/mosaicatcher/analysis/results/sept18_U32/sv_probabilities/HWWKWAFXY_HG03683x01_19s004569-1-1/100000_fixed_norm.selected_j0.1_s0.1/probabilities.Rdata"
#p_link = "/home/hoeps/PhD/projects/huminvs/mosaicatcher/analysis/results/sept22_U32/sv_probabilities/HG00733/100000_fixed_norm.selected_j0.1_s0.1/probabilities.Rdata"
#p_link = "/home/hoeps/Desktop/probabilities.Rdata"
#CN_link = "/home/hoeps/PhD/projects/huminvs/mosaicatcher/tracks2/Mapping_Normalization/CN_track_plots/result/nonred_inversions_n32_CN.txt"
#CN_link = "/home/hoeps/PhD/projects/huminvs/mosaicatcher/tracks/tracks_hufsah_21sept_second/result/00733_CN.txt"
#outdir_raw = "./deleteme/"
#labels_link = NULL
#p_link = '/home/hoeps/PhD/projects/huminvs/mosaicatcher/analysis/results/hack_32/sv_probabilities/HWWKWAFXY_HG03683x01_19s004569-1-1/100000_fixed_norm.selected_j0.1_s0.1/probabilities.Rdata'

p_link = opt$file
labels_link = opt$bed
outdir_raw = opt$outdir
CN_link = opt$cn_map


### switch modules on/off
make_bell_bulk = T
make_table_bulk = T
make_bee_bulk = F

make_bell_sc = F
make_table_sc = F
make_bee_sc = F

suppressMessages(dir.create(outdir_raw))

# is input file specified?
if (is.null(opt$file)){
  print_help(opt_parser)
  stop("Please specify path to probabilities.R", call.=FALSE)
}

#Read the CN_mapability file first
if (is.null(CN_link)){
  print('No mapability/copy number track specified. Skipping that part.')
  CN = NULL
} else {
  print('Mapability/CN file provided. Reading and storing the information')
  CN = read.table(CN_link, stringsAsFactors = F, header=1);
  colnames(CN) = c('chrom','start','end','sample','INV', 'CN', 'mapability', 'n_mappable_bins')
}


# Tell the user if we go for bed mode or single mode. At the occasion, also load it.
if (is.null(labels_link)){
  print('No bed file specified. Examining everything together.')
  labels = NULL
  } else {
  print('Bed file provided. Will use it to stratify results w.r.t. groups')
  labels = read.table(labels_link, stringsAsFactors = FALSE); colnames(labels) = c('chrom','start','end','group')
  
  # Bed files can be weird. For safety, we remove duplicate lines
  labels = unique(labels)
}

# sample name could be interesting
sname = substring(p_link,data.frame(str_locate(p_link, 'HG|NA|GM'))$start,data.frame(str_locate(p_link, 'HG|NA|GM'))$start+6)


# load p to p_grouped
print('Loading probabilities table')
probs_raw = load_and_prep_pdf(p_link)

if (!is.null(CN)){
  CNmerge = as.data.frame(lapply(CN[, c("chrom","start","end","CN","mapability","n_mappable_bins")], as.character))
  CNmerge = as.tbl(CNmerge)
  CNmerge <- CNmerge %>%  mutate(chrom = as.character(chrom),
                                 start = as.numeric(as.character(start)),
                                 end = as.numeric(as.character(end)),
                                 CN = as.numeric(as.character(CN)),
                                 mapability = as.numeric(as.character(mapability)))
  p2 <- full_join(probs_raw, CNmerge, by = c("chrom","start","end"))
  
  # TODO fix
  if (length(p2[is.na(p2$CN),]$CN) > 0){
  p2[is.na(p2$CN),]$CN = 0
  }
  if (length(p2[is.na(p2$mapability),]$mapability) > 0){
  p2[is.na(p2$mapability),]$mapability = 0
  }
}
probs_raw = p2[!is.na(p2$sample),]

len_normalization = as.numeric(as.character(probs_raw$n_mappable_bins))/((probs_raw$end - probs_raw$start)/100.)
probs_raw$expected = probs_raw$expected * len_normalization
#probs_raw$W = probs_raw$W * len_normalization
#probs_raw$C = probs_raw$C * len_normalization

# 
# aaa = (p2[1:10,])
# probs_raw$expected = probs_raw$expected * 

# Adding group information to probs_raw.
if (is.null(labels)){
  # If bed file was not provided, everyone is group 'all'
  probs_raw$group = 'all'
} else {
  # Else, the ones with a group get that one from the bed file
  probs_raw = full_join(probs_raw, labels)
  # ... the remaining ones are called ungrouped.
  probs_raw$group[is.na(probs_raw$group)] = 'ungrouped'
}

# remove inf things [DIRTY SOLUTION! SHOULD BE DONE BETTER! Not needed apparently? Not sure. CHECK]
#probs_raw = probs_raw[probs_raw$logllh != 'Inf',]
#probs_raw = probs_raw[probs_raw$logllh != '-Inf',]
#probs_raw = probs_raw[!(probs_raw$cell %in% unique(probs_raw[probs_raw$logllh == '-Inf',]$cell)),]

###################################
##### OKAAYYY HERE WE GOOOOO ######
#### THIS IS THE MAIN WORKHORSE ###
###################################

group = unique(probs_raw$group)[1] #for quick manual mode
for (group in unique(probs_raw$group)){
  # Talk to human
  print(paste0('Running samples with group ', group))
  
  # Make the outfolder (maybe not necessary?)
  outdir = gsub("\\.:", "_:", paste0(outdir_raw ,group,'/'))
  suppressMessages(dir.create(outdir))
  
  # cut pg down to the desired inversions
  pg = probs_raw[probs_raw$group == group,]
  haps_to_consider = na.omit(unique(pg$haplotype))
  

  if (make_bell_bulk){
    ### [I]a) make dumbbell plot ###
    ####[I] BULK ###
    
    pg_bulk_list = (bulkify_pg(haps_to_consider, pg))
    pg_bulk = data.frame(pg_bulk_list[1]) %>% group_by(start, end, haplotype, class, group)
    pg_bulk_probs = data.frame(pg_bulk_list[2]) %>% group_by(start, end, haplotype, class, group)
    #write.table(pg_bulk, file=paste0('/home/hoeps/Desktop/', 'counts_bulk.txt'), quote = F, row.names = F, col.names = T)
    
    # at least temporarily, I'm operating both with likelihoods and probabilities. Haven't decided yet which one I like more.
    #call_llhs_bulk = (make_condensed_sumlist(haps_to_consider, pg_bulk)) %>% mutate_all(funs(replace_na(.,-1000)))
    call_llhs_bulk = (make_condensed_sumlist(haps_to_consider, pg_bulk)) %>% mutate_all(funs(replace_na(.,-1000)))
    
    #call_probs_bulk = (make_condensed_sumlist_probs(haps_to_consider, pg_bulk_probs)) %>% mutate_all(funs(replace_na(.,-1000)))
    #write.table(mm, file=paste0('/home/hoeps/Desktop/', 'counts_bulk2.txt'), quote = F, row.names = F, col.names = T)
    

    
    g = make_dumbbell(call_llhs_bulk, groupname=group, run_shiny=F)
    p = suppressMessages(ggplotly(g))
  
    savepath = paste0(outdir, 'bellplot_bulk.html')
    suppressMessages(htmlwidgets::saveWidget(as_widget(p), file.path(normalizePath(dirname(savepath)),basename(savepath))))
    ggsave(filename=paste0(outdir, sname, '_', group, '_bellplot_bulk.png'), width=30, height=12, units='cm', device='png')
    ggsave(filename=paste0(outdir, sname, '_', group, '_bellplot_bulk.pdf'), width=30, height=12, units='cm', device='pdf')
    
    #call_probs_bulk[,4:73] = -log(1-(call_probs_bulk[,4:73]))
    #call_probs_bulk[,4:73] = 10**(call_probs_bulk[,4:73])
    #g = make_dumbbell_probs(call_probs_bulk, groupname=group, run_shiny=F)
    #p = suppressMessages(ggplotly(g))
    
    savepath = paste0(outdir, 'bellplot_bulk_prob.html')
    suppressMessages(htmlwidgets::saveWidget(as_widget(p), file.path(normalizePath(dirname(savepath)),basename(savepath))))
    ggsave(filename=paste0(outdir, sname, '_', group, '_bellplot_bulk_prob.png'), width=30, height=12, units='cm', device='png')
    ggsave(filename=paste0(outdir, sname, '_', group, '_bellplot_bulk_prob.pdf'), width=30, height=12, units='cm', device='pdf')
    
  }

  if (make_table_bulk){
    #### [I]b) write table ####
    tab = make_table_finaledition(call_llhs_bulk, group, sname)
    #write.table(tab, file=paste0('/home/hoeps/Desktop/', 'counts_bulk_labels.txt'), quote = F, row.names = F, col.names = T)
    #adding copy number and mapability information to the table
    
    
    t2 = as.data.frame(lapply(tab, as.character))
    
    if (!is.null(CN)){
      CNmerge = as.data.frame(lapply(CN[, c("chrom","start","end","CN","mapability")], as.character))
      tab <- left_join(t2, CNmerge, by = c("chrom","start","end"))
      CN = read.table(CN_link, stringsAsFactors = F, header=1);
    } 
    

    write.table(tab, file=paste0(outdir, 'sv_calls_bulk.txt'), quote = F, row.names = F, col.names = T)
  }
  
  if (make_bee_bulk){
    ## [I]c) make beeswarm plots ##
    print('sup')
    save_beeswarms(pg_bulk %>% group_by(start), call_llhs_bulk, outdir, testrun=F, compositemode=T)
    print('over the hill')
  }

  ###[II] SINGLE CELL ###
  
  
  if (make_bell_sc){
    call_llhs = (make_condensed_sumlist(haps_to_consider, pg))# %>% mutate_all(funs(replace_na(.,-1000)))
    
    #### [II]a) make dumbbell plot ####
    
    
    suppressMessages(source("regenotype_helpers.R")) #for quick manual mode
    
    # create the ggplot plot
    g = make_dumbbell(call_llhs, groupname=group, run_shiny=F)
    g
    # convert it to plotly
    p = suppressMessages(ggplotly(g))
  
    # save. htmlwidgets does not work with relative paths, so we need a little workaround
    # (info and code taken from https://stackoverflow.com/questions/41399795/savewidget-from-htmlwidget-in-r-cannot-save-html-file-in-another-folder)
    savepath = paste0(outdir, 'bellplot.html')
    suppressMessages(htmlwidgets::saveWidget(as_widget(p), file.path(normalizePath(dirname(savepath)),basename(savepath))))
    ggsave(filename=paste0(outdir, sname,'_', group,'_bellplot.png'), width=30, height=12, units='cm', device='png')
    ggsave(filename=paste0(outdir, sname,'_', group,'_bellplot.pdf'), width=30, height=12, units='cm', device='pdf')
    }
    
  if (make_table_sc){
    #### [II]b) write table ####
    #tab = make_table(call_llhs, group, sname)
    tab = make_table_finaledition(call_llhs, group, sname)
    
    t2 = as.data.frame(lapply(tab, as.character))
    CNmerge = as.data.frame(lapply(CN[, c("chrom","start","end","CN","mapability")], as.character))
    tab2 <- left_join(t2, CNmerge, by = c("chrom","start","end"))
    
    
    #adding copy number and mapability information to the table
    write.table(tab2, file=paste0(outdir, 'sv_calls.txt'), quote = F, row.names = F, col.names = T)
  }
  
  if (make_bee_sc){
    #### [II]c) save beewarm plots ###
    save_beeswarms(pg, call_llhs, outdir, testrun=F)
    print(paste0('Group ', group,' done.'))
  }
}
print('### ALL DONE. Happy discoveries. ###')

# Does not seem to work.
options(warn = oldw)



# Manual links for debugging. 
#p_link = '/home/hoeps/PhD/projects/huminvs/mosaicatcher/sv_probabilities/HG00733/100000_fixed_norm.selected_j0.1_s0.5/probabilities.Rdata'
# p_link ='../results/audano_calls/probabilities.Rdata'
# outdir = '/home/hoeps/PhD/projects/huminvs/mosaicatcher/analysis/results/31st_may/'

