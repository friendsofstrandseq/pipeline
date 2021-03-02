# Whoeps, 07th Jan 2021
# Making a large overview over the results from the arbigent folder. 
# I'm giving myself 2h to make this nice today. 

# Input: callmatrix from clean_genotype.R
# Input: a csv from david from which to extract samplenames
# Output: a matrix with added entries:
#   Filter - Pass, NoReadsPass, MendelFail, FalsePositive, lowconf, AlwaysComplex
#   Mapability - percent?
#   nhom, nhet, nref, nnoreads, ncomplex
#   mendel 1/0

# Load libraries
library(stringr)
library(dplyr)
library(pheatmap)
library(matrixStats)
library(reshape2)


### FUNCTIONS ###


extract_used_samples <- function(csvlink){
  # This is for processing of one of David's typical outputs. #
  # Load david's genotypes
  
  dgt = read.table(david_list, header=T, sep=',', stringsAsFactors = F)
  
  # define entries to keep
  cols_to_keep = c('genoT')
  
  # filter
  dgtf = dgt %>%   select(matches(paste(cols_to_keep, collapse="|")))
  
  # rename columns
  colnames(dgtf) = str_replace(colnames(dgtf),'genoT_','')
  
  return(colnames(dgtf))
  
}

cut_cm_to_samples <- function(cm_f, used_samples_uniq_f, n_overjump, include_children){
  
  # Kick children out if not wanted (hehe)
  if (include_children){
    children_ids = c('HG00514','HG00733','NA19240')
  } else {
    children_ids = NULL
  }
  # Change callmatrix samplenames
  
  # Go in reverse order throu sample columns
  for (i in length(colnames(cm_f)):n_overjump){
    
    # Rename HG002 if you encounter it
    if (colnames(cm_f)[i]=='HG002'){
      colnames(cm_f)[i] = 'NA24385'
      next
    } else if (colnames(cm_f)[i] %in% children_ids){
      next
    }
    
    # Use the samplenames from David's table to adjust ours. 
    newlabel = used_samples_uniq_f[grepl(substr(colnames(cm_f)[i],3,8),used_samples_uniq_f)]
    print(newlabel)
    # If the sample wasn't found in David's table, we drop it. 
    if (length(newlabel)>0){
      colnames(cm_f)[i] = newlabel
    } else {
      print(paste0('Dropped sample ',colnames(cm)[i]))
      cm_f[,colnames(cm_f)[i]] = NULL
    }
    # Reset newlabel
    newlabel = NULL
  }
  
  return(cm_f)
}

count_homhetrefetc <- function(cm_f, n_samples){
  cm_f$nhom = rowSums2(cm_f=='1|1') + rowSums2(cm_f=='1|1_lowconf')
  cm_f$nhet = rowSums2(cm_f=='1|0') + rowSums2(cm_f=='0|1') + rowSums2(cm_f=='1|0_lowconf') + rowSums2(cm_f=='0|1_lowconf')
  cm_f$nref = rowSums2(cm_f=='0|0') + rowSums2(cm_f== '0|0_lowconf')
  cm_f$nnoreads = rowSums2(cm_f=='noreads')
  cm_f$inv_dup_het <- rowSums2(cm_f == "1011") + rowSums2(cm_f=='1011_lowconf') + rowSums2(cm_f == "1110") + rowSums2(cm_f=='1110_lowconf')
  cm_f$inv_dup_hom <- rowSums2(cm_f == "1111") + rowSums2(cm_f=='1111_lowconf')+ rowSums2(cm_f == "2200") + rowSums2(cm_f=='2200_lowconf')+ rowSums2(cm_f == "0022") + rowSums2(cm_f=='0022_lowconf')
  cm_f$ncomplex = n_samples - (cm_f$nhom + cm_f$nhet + cm_f$nref + cm_f$nnoreads+ cm_f$inv_dup_hom + cm_f$inv_dup_het)
  return(cm_f)
}

apply_filter_new <- function(cm_f, samples){
  
  # Judgement day: filter
  cm_f$gt_events = cm_f$nhom + cm_f$nhet #+ cm_f$ncomplex 
  cm_f$inv_dups <- cm_f$inv_dup_het + cm_f$inv_dup_hom
  cm_f$verdict = 'UNK'
  bincutoff = 5
  
  cols_to_add = c('lowconf','noreads','FP','alwayscomplex','MISO','Mendelfail', 'INVDUP')
  
  # Our criteria: 

  # Lowconf if below bincutoff
  cm_f$lowconf = cm_f$valid_bins <= bincutoff
  
  # noreads if noreads
  cm_f$noreads = cm_f$nnoreads == length(samples)
  
  # FP: all samples are REF
  cm_f$FP = cm_f$nref == length(samples)
  
  # alwayscomplex: all samples are complex (complex includes invdup)
  cm_f$alwayscomplex= cm_f$ncomplex == length(samples)
  
  # MISO: >80% of GTs are hom inversions
  cm_f$MISO = cm_f$nhom >= 0.8*length(samples)
  
  # Mendelfail: At least one trio mendelfails
  cm_f$Mendelfail = cm_f$mendelfails > 0
  
  # INVDUP: any INVDUP observed?
  cm_f$INVDUP = cm_f$inv_dups>0
  
  # Now, merge all these criteria into one column that lists all applying ones. 
  cm_f_i = cm_f[,cols_to_add]
  
  cm_f$verdict = apply(cm_f_i, 1, function(x) paste(names(cm_f_i)[x == T], collapse = "-"))
  cm_f[cm_f$verdict == '',]$verdict = 'pass'
  
  # Clean up
  cm_f[,cols_to_add] = NULL
  cm_f$gt_events = NULL
  return(cm_f)
}


apply_filter <- function(cm_f, used_samples_uniq){
  # Judgement hour: filter
  cm_f$gt_events = cm_f$nhom + cm_f$nhet + cm_f$ncomplex 
  
  # Criteria hour
  # If we have < 5 bins, we do not attempt to filter
  
  # Failure modes
  # FalsePositive: no nonref events confirmed
  # AlwaysComplex: no noncomplex events confirmed
  
  # Explicit pass modes
  # Undercalled: CR is below 0.5, and we observe more events in GT than in the discovery list.
  # Pass: Events match perfectly 
  # NoReadsPass: Events match perfectly, difference can be explained by nnoreads (?)
  # PassMulticallers 
  
  # Judgement hour
  cm_f$verdict = 'UNK'
  bincutoff = 5
  cm_f[(cm_f$valid_bins <= bincutoff),]$verdict = 'lowconf'
  cm_f[(cm_f$gt_events == 0) & (cm_f$valid_bins > bincutoff),]$verdict = 'NoEvents'
  cm_f[(cm_f$mendelfails > 0) & (cm_f$valid_bins > bincutoff),]$verdict = 'MendelFail'
  cm_f[((cm_f$nhom + cm_f$nhet + cm_f$nref) == 0) & (cm_f$valid_bins > bincutoff),]$verdict = 'AlwaysComplex'
  
  cm_f[(cm_f$nnoreads == length(used_samples_uniq)),]$verdict = 'NoReadsPass'
  #cm_f[(cm_f$gt_events == cm_f$eventsclaimed),]$verdict = 'Pass'
  cm_f[cm_f$verdict == 'UNK',]$verdict = 'Pass'
  #cm_f[cm_f$CALLERSET_LIST %in% c('PAV,SSEQAUTO','PAV,SSEQAUTO,BIONANO','SSEQAUTO,BIONANO'),]$verdict = 'Pass_multicallers'
  return(cm_f)
}

#' Make a vector to check inheritance plaisibilities
#' @return child_expect, a vector describing excpected child genotypes given the parents.
#' @author Wolfram Hoeps
make_child_expect_vector <- function(){
  # ok who cares let's do it the hard way.
  #parents = paste(p1, p2)
  child_expect <- vector(mode="list")
  child_expect[['0|0 0|0']] = c('0|0')
  child_expect[['0|0 0|1']] = c('0|0','0|1')
  child_expect[['0|0 1|1']] = c('0|1')
  child_expect[['0|1 0|0']] = c('0|0','0|1')
  child_expect[['0|1 0|1']] = c('0|0','0|1','1|1')
  child_expect[['0|1 1|1']] = c('0|1','1|1')
  child_expect[['1|1 0|0']] = c('0|1')
  child_expect[['1|1 0|1']] = c('0|1','1|1')
  child_expect[['1|1 1|1']] = c('1|1')
  return (child_expect)
}

#' Test a gt for validity
#' 
test_mendel <- function(ce, gt_parent1, gt_parent2, gt_child){
  gt_parent1 = substr(gt_parent1,1,3)
  gt_parent2 = substr(gt_parent2,1,3)
  gt_child = substr(gt_child, 1,3)
  gt_parent1 = gsub('1\\|0','0\\|1', gt_parent1)
  gt_parent2 = gsub('1\\|0','0\\|1', gt_parent2)
  gt_child   = gsub('1\\|0','0\\|1', gt_child)
  
  valid_gts = c('0|0','0|1','1|1')
  c1 = gt_parent1 %in% valid_gts
  c2 = gt_parent2 %in% valid_gts
  c3 = gt_child %in% valid_gts
  #print(c1)
  #print(c2)
  #print(c3)
  #return(gt_parent2)
  if (c1){
    if (c2){
      if (c3){
        valid = gt_child %in% ce[[paste(gt_parent1, gt_parent2)]]
        #print(valid)
        return(valid)   
      }
    }
  }
  return('CPX')
}


add_mendelfails <- function(cm_f){
  ce = make_child_expect_vector()
  callmatrix = cm_f
  callmatrix$mendel1 = 'UNK'
  callmatrix$mendel2 = 'UNK'
  callmatrix$mendel3 = 'UNK'
  
  for (row in 1:nrow(callmatrix)){
    #callmatrix[row,]$mendel = 
    callmatrix[row,]$mendel1 = (test_mendel(ce, callmatrix[row,]$NA19238, callmatrix[row,]$NA19239,callmatrix[row,]$NA19240 ))
    callmatrix[row,]$mendel2 = (test_mendel(ce, callmatrix[row,]$HG00512, callmatrix[row,]$HG00513,callmatrix[row,]$HG00514 ))
    callmatrix[row,]$mendel3 = (test_mendel(ce, callmatrix[row,]$HG00731, callmatrix[row,]$HG00732,callmatrix[row,]$HG00733 ))
    #callmatrix[row,]$mendel4 = as.logical(test_mendel(ce, callmatrix[row,]$GM19650, callmatrix[row,]$HG00864,callmatrix[row,]$HG03371 ))
    
  }
  
  callmatrix$mendelfails = rowSums(callmatrix[,c('mendel1','mendel2','mendel3')] == 'FALSE')
  return(callmatrix)
}
