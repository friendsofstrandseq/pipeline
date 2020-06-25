# Whoeps, 31 May 2020
# ### Preface: This program has probably still 500 errors. 
# ### if something is not running, it's most likely our fault
# ### and not yours. So don't waste your time and ask! :) 
# Input: - probabilities.R
#        - bed file with groups
#        - outdir
#        
# Output: - table with classifications
#         - dumbbellplot
#         - beeswarmplots

# supposed to suppress warnings, but not working
oldw <- getOption("warn")
options(warn = -1)

#!/usr/bin/env Rscript
suppressMessages(library("optparse"))
suppressMessages(source("probability_helpers.R"))
suppressMessages(source("regenotype_helpers.R"))

# INPUT INSTRUCTIONS
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="probabilities.R, produced by mosaicatcher", metavar="character"),
  make_option(c("-b", "--bed"), type="character", default=NULL, 
              help="a bed file specifiying labels/groups for the segments.", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default="./output/", 
              help="output dir name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# is input file specified?
if (is.null(opt$file)){
  print_help(opt_parser)
  stop("Please specify path to probabilities.R", call.=FALSE)
}

# Ok lets go! Calm the minds of impatient humans first of all.
print('Processing and summarizing information, making plots')
print('This can take a few minutes.')

### This block is for debug mode. ###
#p_link = '/home/hoeps/PhD/projects/huminvs/mosaicatcher/sv_probabilities/HG00733/100000_fixed_norm.selected_j0.1_s0.5/probabilities.Rdata'
#p_link ='../../results/HG733_postcorr2/probabilities.Rdata'
#outdir_raw = 'aud2/'
#labels_link = '/home/hoeps/PhD/projects/huminvs/mosaicatcher/analysis/results/HG733_postcorr/naming_HG733.bed'
#####################################

p_link = opt$file
xlabels_link = opt$bed
outdir_raw = opt$outdir



suppressMessages(dir.create(outdir_raw))


# Tell the user if we go for bed mode or single mode. At the occasion, also load it.
if (is.null(labels)){
  print('No bed file specified. Examining everything together.')
  } else {
  print('Bed file provided. Will use it to stratify results w.r.t. groups')
  labels = read.table(labels_link, stringsAsFactors = FALSE); colnames(labels) = c('chrom','start','end','group')
}



# load p to p_grouped
pg_a = load_and_prep_pdf(p_link)
if (is.null(labels)){
  # If bed file was not provided, everyone is group 'all'
  pg_a$group = 'all'
} else {
  # Else, the ones with a group get that one from the bed file
  pg_a = full_join(pg_a, labels)
  # ... the remaining ones are called ungrouped.
  pg_a$group[is.na(pg_a$group)] = 'ungrouped'
}

# if bed file contained inversions not present in probabilities, we need this.
pg_a = na.omit(pg_a)
###################################
##### OKAAYYY HERE WE GOOOOO ######
#### THIS IS THE MAIN WORKHORSE ###
###################################

#go for each group separately
# group = unique(pg_a$group)[1] #for quick manual mode
for (group in unique(pg_a$group)){
  
  # Talk to human
  print(paste0('Running samples with group ', group))
  
  # Make the outfolder (maybe not necessary?)
  outdir = paste0(outdir_raw ,group,'/')
  suppressMessages(dir.create(outdir))
  
  # cut pg down to the desired inversions
  pg = pg_a[pg_a$group == group,]
  
  haps_to_consider = na.omit(unique(pg$haplotype))
  
  # and here is in essence already our result. 
  call_llhs = make_condensed_sumlist(haps_to_consider, pg)
  
  
  #### a) make dumbbell plot ####
  
  #suppressMessages(source("standalone_functions.R")) #for quick manual mode
  
  # create the ggplot plot
  g = make_dumbbell(call_llhs, run_shiny=F)
  g
  # convert it to plotly
  p = suppressMessages(ggplotly(g))

  # save. htmlwidgets does not work with relative paths, so we need a little workaround
  # (info and code taken from https://stackoverflow.com/questions/41399795/savewidget-from-htmlwidget-in-r-cannot-save-html-file-in-another-folder)
  savepath = paste0(outdir, 'bellplot.html')
  suppressMessages(htmlwidgets::saveWidget(as_widget(p), file.path(normalizePath(dirname(savepath)),basename(savepath))))
  ggsave(filename=paste0(outdir, 'bellplot.png'), width=30, height=12, units='cm', device='png')
  
  #### b) write table ####
  tab = make_table(call_llhs)
  write.table(tab, file=paste0(outdir, 'sv_calls.txt'), quote = F, row.names = F, col.names = T)
  
  #### c) save beewarm plots ###
  save_beeswarms(pg, call_llhs, outdir, testrun=F)
  print(paste0('Group ', group,' done.'))
}
print('### ALL DONE. Happy discoveries. ###')

# Does not seem to work.
options(warn = oldw)



# Manual links for debugging. 
#p_link = '/home/hoeps/PhD/projects/huminvs/mosaicatcher/sv_probabilities/HG00733/100000_fixed_norm.selected_j0.1_s0.5/probabilities.Rdata'
# p_link ='../results/audano_calls/probabilities.Rdata'
# outdir = '/home/hoeps/PhD/projects/huminvs/mosaicatcher/analysis/results/31st_may/'
