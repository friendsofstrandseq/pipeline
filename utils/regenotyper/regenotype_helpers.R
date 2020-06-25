# Whoeps, 31th May 2020
# Admittedly the name 'standalone' is not wisely chosen. Anyway these are functions needed for the 
# three outputs of standalone. R

#' Function A: make a dumbbell plot
#'
#' @param call_llhs
#' @param 
#' @author Wolfram Hoeps
#' @export
make_dumbbell <- function(segs_llhs_f, run_shiny=F){
  library(ggplot2)
  library(plotly)
  library(shiny)
  library(RColorBrewer)
  library(scales)
  # take the segs_llhs (which are call_llhs), and add important information:
  # maxval, secval, len, maxname, interesting_minval and sum_hetinv_max
  segs_llhs = attach_max_sec_name_to_call_llhs(segs_llhs_f)
  
  # I think len is no longer needed but I could be wrong.
  segs_llhs$len = segs_llhs$end - segs_llhs$start
  segs_llhs$interesting_minval = sign(segs_llhs$interesting_minval) * log10(abs(segs_llhs$interesting_minval)+1)
  segs_llhs$maxval = sign(segs_llhs$maxval) *log10(abs(segs_llhs$maxval)+1)
  segs_llhs$sum_hetinv_max = sign(segs_llhs$sum_hetinv_max) *log10(abs(segs_llhs$sum_hetinv_max)+1)
  segs_llhs$sum_0101 = sign(segs_llhs$sum_0101) *log10(abs(segs_llhs$sum_0101)+1)
  segs_llhs$maxname_short = substr(segs_llhs$maxname,5,9)
  ### PLOT ###
  g <- ggplot(segs_llhs, aes(x=(len/1000))) +scale_fill_brewer(palette = "Spectral")
  g = g  +    
    # Add a few grid lines first of all
    geom_hline(yintercept=log1p(0), alpha=1) +
    geom_hline(yintercept=2, alpha=0.2) +
    geom_hline(yintercept=-2, alpha=0.2) +
    
    # Add the vertical lines that go from the lowest to the highest value for each segment
    geom_segment(aes(xend=(len/1000), y=interesting_minval, yend=maxval), 
                 size=0.5,  colour="grey70") + 
    
    # Add the dots for HOM, HET and HIGHEST. They will all lie on the vertial lines we just created
    geom_point(aes(y=maxval, colour='score OTHER', group=segs_llhs$start, text=maxname_short), 
               alpha=1, size=2.5) +#, col=brewer.pal(n=3, name='Set1')[1]) +
    geom_point(aes(y=sum_hetinv_max, colour='score HET',group=segs_llhs$start),
               alpha=1, size=2.5) +#, col=brewer.pal(n=3, name='Set1')[3]) +
    geom_point(aes(y=sum_0101, colour='score HOM',group=segs_llhs$start),
               alpha=1, size=2.5) +#, col=brewer.pal(n=3, name='Set1')[2]) +

    # Change colors. Entries are sorted alphabetically, so this is a bit a mess. 
    # alphabet: [het, hom, other], hue_pal: red, green, blue
    # desired order: het - green, other - red, hom - blue
    # therefore: c(2,3,1)
    scale_color_manual(values=hue_pal()(3)[c(2,3,1)]) +
    
    # Add labels
    labs(title="Inversion re-genotyping", 
         subtitle="") +
    xlab('SV length [kb]') +
    ylab('LLR over REF') +
    
    # Make a nice x axis with log and logticks
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x)
      #limits=c(1,5000)
      #labels = scales::trans_format("log10", scales::x)
    ) +
    # scale_y_log10(
    #   breaks = scales::trans_breaks("log10", function(y) 10^(y+1))
    #   #limits=c(1,5000)
    #   #labels = scales::trans_format("log10", scales::x)
    # ) +
    annotation_logticks(sides = 'b',colour = 'gray20' ) +
    scale_y_continuous(breaks = (seq(-4,4,by=1))) +
    # Remove title for all legends
    #theme(legend.title=element_blank()) +
    
    # Add text for other
    geom_text(aes(y= maxval, label = maxname_short),
              size = 2,
              hjust = "left") +
    
    # EEEEEHm, not sure. 
    guides(color = guide_legend(reverse = TRUE))
  
  # Okay, we created 'g' now. We CAN pass it to shiny if we want to
  if (run_shiny){
    # embed in plotly #
    ui <- fluidPage(
      plotlyOutput("distPlot")
    )
    server <- function(input, output) {
      output$distPlot <- renderPlotly({
        g
      })
    }
    # Run a shiny app
    shinyApp(ui = ui, server = server)
  }
  
  return(g)
}

#' Function C: save a beeswarm plot for each segment
#'
#' @param call_llhs
#' @param 
#' @author Wolfram Hoeps
#' @export
save_beeswarms <- function(pg_f, call_llhs_f, outdir, testrun=F){
  suppressMessages(dir.create(outdir))
  suppressMessages(dir.create(paste0(outdir, 'beeplots')))
  suppressMessages(dir.create(paste0(outdir, 'beeplots/0101')))
  suppressMessages(dir.create(paste0(outdir, 'beeplots/0110')))
  suppressMessages(dir.create(paste0(outdir, 'beeplots/1001')))
  suppressMessages(dir.create(paste0(outdir, 'beeplots/1010')))
  suppressMessages(dir.create(paste0(outdir, 'beeplots/else')))
  suppressMessages(dir.create(paste0(outdir, 'beeplots/reject')))
  
  only_main_haps=F
  
  # Iterate over all segments
  if (testrun==T){
    print('Running beewarm in testmode. Only first 10 segments considered.')
    n_pics = 3
  } else {
    n_pics = dim(call_llhs)[1]
  }
  for (n in 1:n_pics){
    
    # Extract the information of for this segment
    pgi = suppressMessages(select_group(pg_f, n))
    # Return a list with the most likely classifications
    topnames_hap_presort = get_top_scoring_haplotypes_standalone(call_llhs_f, pgi, 70)
    # ...and extract the absolute winner FIRST
    most_likely_pred = topnames_hap_presort[1]
    # ... and now sort it
    topnames_hap = c(c('0101','1001','0110'), topnames_hap_presort[!topnames_hap_presort %in% c('0101', '1001','0110')])
    # ... and return the 30 top ones
    topnames_hap = topnames_hap[1:30]
    ## In order to make the 'winner cell' larger, we need to do some magic here. ##
    # pgi2 are all data
    pgi2 = pgi[pgi$haplotype %in% topnames_hap,]
    # pgi3 are SV llhs of only the winning SV per cell.
    pgi3 = pgi %>% group_by(cell) %>% top_n(1, logllh)
    pgi3 = pgi3[pgi3$haplotype %in% topnames_hap,]
    
    ## magic end ##
    
    #### PLOT ####
    g <- ggplot(data=pgi2, aes(x=haplotype, y=logllh)) +  
      # Plot first the fat dots
      geom_beeswarm(data=pgi3,
                    size=0.5,
                    cex=0.2,
                    aes(color=class)) +
      # ... then add the slim ones
      geom_beeswarm(data=pgi2,
                    size=0.2,
                    cex=0.1,
                    alpha=0.5,
                    aes(color=class)
      ) +
      # Add median bars
      stat_summary(fun.y=median, fun.ymin=median, fun.ymax=median, 
                   geom="crossbar", width=0.7, size=0.1) +
      
      # Sort entries by our topnames list
      scale_x_discrete(limits = topnames_hap) +
      theme(axis.text.x = element_text(angle = 90)) +

      # Add title
      labs(title=paste(pgi2$chrom, ':',pgi2$start, '-',pgi2$end, ' ', ' LEN=',pgi2$len/1000, '  ', ' n=', n, sep=''), subtitle="") +
      
      # Dotted line for 0
      geom_hline(yintercept = 0, linetype="dashed", 
                 color = "black", size=0.3) 
    
    # Save to the outbeedir
    if (most_likely_pred %in% c('1010','0101','0110','1001')){
      ggsave(filename=paste(outdir, 'beeplots/',most_likely_pred,'/', pgi2$chrom[1], ':',pgi2$start[1], '-',pgi2$end[1], '_',round(pgi2$len[1]/1000, 3),'k.png',sep=''), width=30, height=12, units='cm', device='png')
    } else {
      if (test_for_invrejection(call_llhs_f, pgi2$chrom[1], pgi2$start[1])){
        # If TRUE was returned, it means the inversion can be rejected,
        # i.e. HOM and HET are both less likely than reference. 
        ggsave(filename=paste(outdir, 'beeplots/','reject','/', pgi2$chrom[1], ':',pgi2$start[1], '-',pgi2$end[1], '_',round(pgi2$len[1]/1000, 3),'k.png',sep=''), width=30, height=12, units='cm', device='png')
      } else {
        # These are all segments that have no strongest HOM or HET fighter, but also HOM and HET are not categorically out of question yet. This will be a mixture of light green, red and grey.
        ggsave(filename=paste(outdir, 'beeplots/','else','/', pgi2$chrom[1], ':',pgi2$start[1], '-',pgi2$end[1], '_',round(pgi2$len[1]/1000, 3),'k.png',sep=''), width=30, height=12, units='cm', device='png')
      }
    }
    print(paste0('Saved number ', n))
  }
}




#' Function B: return an overview table
#'
#' @param call_llhs
#' @author Wolfram Hoeps
#' @export
make_table <- function(call_llhs_f){
  
  # take the segs_llhs (which are call_llhs), and add important information:
  # maxval, secval, len, maxname, interesting_minval and sum_hetinv_max
  call_llhs_f2 = attach_max_sec_name_to_call_llhs(call_llhs_f)
  
  # also the length. it's a bit a messy topic.
  call_llhs_f2$len = (call_llhs_f2$end - call_llhs_f2$start)/1000

  # simply extract the interesting information
  cols_of_interest = c('chrom', 'start', 'end', 'len', 'maxname','secname','maxval', 'secval')
  overview = call_llhs_f2[,cols_of_interest]
  
  # round some numbers
  overview$maxval = round(overview$maxval, 1)
  overview$secval = round(overview$secval, 1)
  
  # sort
  overview = overview[with(overview, order(maxname, -maxval)), ]

  return(overview)
}

#' Little helperfunction to decide if an inversion can be automatically rejected. This is the case if
#' the LLHs for HOM and both HETs are below reference. It wouldn't have been necessary to write its own
#' function for this, but well here it is now.
#'
#' @param cl call_llhs
#' @param chrom the chromosome of the inversion of interest
#' @param start the start position of the inversion of interest
#' @author Wolfram Hoeps
#' @export
test_for_invrejection <- function(cl, chrom, start){
  # Read the values for inv_hom, inv_h1, inv_h2
  vals_to_check = cl[(cl$chrom==chrom) & (cl$start==start),][c('sum_0101', 'sum_0110', 'sum_1001')]
  # Are all of them < 0? Then we can reject.
  can_be_rejected = all(vals_to_check < 0)
  # True: reject. False: keep
  return(can_be_rejected)
}