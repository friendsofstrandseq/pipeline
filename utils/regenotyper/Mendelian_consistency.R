# Hufsah Ashraf, 4-2-2021
#this function takes per sample inversion genotypes and checks if they are Mendelian Consistent for each of the 3 trios

mendel_cons<-function(genotypes){
library(stringr)
library(stringi)
library(dplyr)
mendel<- data.frame(genotypes$chrom, genotypes$start, genotypes$end, genotypes$HG00731, genotypes$HG00732, genotypes$HG00733,
                    genotypes$NA19239, genotypes$NA19238, genotypes$NA19240, genotypes$HG00512, genotypes$HG00513, genotypes$HG00514)
trio_1_genots<-data.frame(mendel$genotypes.chrom, mendel$genotypes.start, mendel$genotypes.end, mendel$genotypes.HG00512, mendel$genotypes.HG00513,
                          mendel$genotypes.HG00514)
trio_1_genots$sample<-'HG00514'
colnames(trio_1_genots)<-c('chrom','start', 'end', 'father', 'mother', 'child', 'sample')
trio_2_genots<-data.frame(mendel$genotypes.chrom, mendel$genotypes.start, mendel$genotypes.end, mendel$genotypes.HG00731, mendel$genotypes.HG00732,
                          mendel$genotypes.HG00733)
trio_2_genots$sample<-'HG00733'
colnames(trio_2_genots)<-c('chrom','start', 'end', 'father', 'mother', 'child', 'sample')
trio_3_genots<-data.frame(mendel$genotypes.chrom, mendel$genotypes.start, mendel$genotypes.end, mendel$genotypes.NA19238, mendel$genotypes.NA19239,
                          mendel$genotypes.NA19240)
trio_3_genots$sample<-'NA19240'
colnames(trio_3_genots)<-c('chrom','start', 'end', 'father', 'mother', 'child', 'sample')
trio1_2_genots<-rbind(trio_1_genots,trio_2_genots)
all_trios_genots<-rbind(trio1_2_genots,trio_3_genots)
#make genotypes simpler to make  comparison easier
mendel<- data.frame(lapply(mendel, function(x) { gsub("_lowconf", "", x)  }))
mendel<- data.frame(lapply(mendel, function(x) { gsub(".\\|\\.", "noreads", x)  }))
mendel<- data.frame(lapply(mendel, function(x) { str_replace_all(x, ".\\/.", "noreads")  }))
mendel$HG00733_cons<-NA
mendel$NA19240_cons<-NA
mendel$HG00514_cons<-NA

#encode 2 letter genotypes To 4 letter codes
mendel<- data.frame(lapply(mendel, function(x) {str_replace_all(x, '0\\|0', "1010") }))
mendel<- data.frame(lapply(mendel, function(x) {str_replace_all(x, '1\\|0', "0110")  }))
mendel<- data.frame(lapply(mendel, function(x) {str_replace_all(x, '0\\|1', "1001") }))
mendel<- data.frame(lapply(mendel, function(x) {str_replace_all(x, '1\\|1', "0101")}))

#exclude 'noreads' regions from testing
mendel$HG00733_cons<-ifelse((as.character(mendel$genotypes.HG00731)=='noreads'|
                               as.character(mendel$genotypes.HG00732)=='noreads'| 
                               as.character(mendel$genotypes.HG00733)=='noreads'),'noreads', NA)


mendel$NA19240_cons<-ifelse((as.character(mendel$genotypes.NA19239)=='noreads'|
                               as.character(mendel$genotypes.NA19238)=='noreads'|
                               as.character(mendel$genotypes.NA19240)=='noreads'), 'noreads', NA)


mendel$HG00514_cons<-ifelse((as.character(mendel$genotypes.HG00513)=='noreads'|
                               as.character(mendel$genotypes.HG00512)=='noreads'| 
                               as.character(mendel$genotypes.HG00514)=='noreads'),'noreads', NA)


#start testing
#HG00514
complex<- data.frame(mendel[is.na(mendel$HG00514_cons),]) 
complex$HG00514_H1_M<-as.character(stri_sub_all(complex[,'genotypes.HG00513'], 1, 2))
complex$HG00514_H2_M<-as.character(stri_sub_all(complex[,'genotypes.HG00513'], 3, 4))
complex$HG00514_H1_F<-as.character(stri_sub_all(complex[,'genotypes.HG00512'], 1, 2))
complex$HG00514_H2_F<-as.character(stri_sub_all(complex[,'genotypes.HG00512'], 3, 4))
complex$HG00514_H1_C<-as.character(stri_sub_all(complex[,'genotypes.HG00514'], 1, 2))
complex$HG00514_H2_C<-as.character(stri_sub_all(complex[,'genotypes.HG00514'], 3, 4))
mendel[is.na(mendel$HG00514_cons),]$HG00514_cons<-ifelse(((complex$HG00514_H1_C==complex$HG00514_H1_M | complex$HG00514_H1_C==complex$HG00514_H2_M) &
                                                            (complex$HG00514_H2_C==complex$HG00514_H1_F | complex$HG00514_H2_C==complex$HG00514_H2_F))|
                                                           ((complex$HG00514_H1_C==complex$HG00514_H1_F | complex$HG00514_H1_C==complex$HG00514_H2_F)&
                                                              (complex$HG00514_H2_C==complex$HG00514_H1_M | complex$HG00514_H2_C==complex$HG00514_H2_M)),
                                                         'YES', 'NO')
#HG00733
complex<- data.frame(mendel[is.na(mendel$HG00733_cons),]) 
complex$HG00733_H1_M<-as.character(stri_sub_all(complex[,'genotypes.HG00732'], 1, 2))
complex$HG00733_H2_M<-as.character(stri_sub_all(complex[,'genotypes.HG00732'], 3, 4))
complex$HG00733_H1_F<-as.character(stri_sub_all(complex[,'genotypes.HG00731'], 1, 2))
complex$HG00733_H2_F<-as.character(stri_sub_all(complex[,'genotypes.HG00731'], 3, 4))
complex$HG00733_H1_C<-as.character(stri_sub_all(complex[,'genotypes.HG00733'], 1, 2))
complex$HG00733_H2_C<-as.character(stri_sub_all(complex[,'genotypes.HG00733'], 3, 4))
mendel[is.na(mendel$HG00733_cons),]$HG00733_cons<-ifelse(((complex$HG00733_H1_C==complex$HG00733_H1_M | complex$HG00733_H1_C==complex$HG00733_H2_M) &
                                                            (complex$HG00733_H2_C==complex$HG00733_H1_F | complex$HG00733_H2_C==complex$HG00733_H2_F))|
                                                           ((complex$HG00733_H1_C==complex$HG00733_H1_F | complex$HG00733_H1_C==complex$HG00733_H2_F)&
                                                              (complex$HG00733_H2_C==complex$HG00733_H1_M | complex$HG00733_H2_C==complex$HG00733_H2_M)),
                                                         'YES', 'NO')
#NA19240
complex<- data.frame(mendel[is.na(mendel$NA19240_cons),]) 
complex$NA19240_H1_M<-as.character(stri_sub_all(complex[,'genotypes.NA19239'], 1, 2))
complex$NA19240_H2_M<-as.character(stri_sub_all(complex[,'genotypes.NA19239'], 3, 4))
complex$NA19240_H1_F<-as.character(stri_sub_all(complex[,'genotypes.NA19238'], 1, 2))
complex$NA19240_H2_F<-as.character(stri_sub_all(complex[,'genotypes.NA19238'], 3, 4))
complex$NA19240_H1_C<-as.character(stri_sub_all(complex[,'genotypes.NA19240'], 1, 2))
complex$NA19240_H2_C<-as.character(stri_sub_all(complex[,'genotypes.NA19240'], 3, 4))
mendel[is.na(mendel$NA19240_cons),]$NA19240_cons<-ifelse(((complex$NA19240_H1_C==complex$NA19240_H1_M | complex$NA19240_H1_C==complex$NA19240_H2_M) &
                                                            (complex$NA19240_H2_C==complex$NA19240_H1_F | complex$NA19240_H2_C==complex$NA19240_H2_F))|
                                                           ((complex$NA19240_H1_C==complex$NA19240_H1_F | complex$NA19240_H1_C==complex$NA19240_H2_F)&
                                                              (complex$NA19240_H2_C==complex$NA19240_H1_M | complex$NA19240_H2_C==complex$NA19240_H2_M)),
                                                         'YES', 'NO')

#verdict for each inversion based on Mendelian consistency check
mendel_genotypes<-data.frame(mendel$genotypes.chrom, mendel$genotypes.start,mendel$genotypes.end, mendel$HG00733_cons, mendel$HG00514_cons,mendel$NA19240_cons)
for (m in 1:length(mendel_genotypes$mendel.genotypes.chrom)){
  #how many trios passed the test
  mendel_genotypes[m,'mendel_cons']<-ifelse(mendel_genotypes[m,'mendel.HG00733_cons']=='noreads'| mendel_genotypes[m,'mendel.HG00514_cons']=='noreads'|
                                              mendel_genotypes[m,'mendel.HG00514_cons']=='noreads' , NA, length(which(mendel_genotypes[m,]=='YES')))
  #whether the inversion passed the mendelian consistency check or not(atleast 2 trios should pass), '0' for Pass and '1' for Fail
  mendel_genotypes[m,'mendel_fail']<-ifelse(mendel_genotypes[m,'mendel_cons']>=2, 0, 1)
}
colnames(mendel_genotypes)<-c('chrom','start','end', 'HG00733_cons','HG00514_cons', 'NA19240_cons','Mendel_cons_trios', 'Mendel_fail')
mendel_genotypes<-data.frame(mendel_genotypes$chrom, mendel_genotypes$start, mendel_genotypes$end, mendel_genotypes$Mendel_cons_trios, mendel_genotypes$Mendel_fail)
colnames(mendel_genotypes)<-c('chrom','start','end','Mendel_cons_trios', 'Mendel_fail')
mendel_genotypes$start<-as.integer(as.character(mendel_genotypes$start))
mendel_genotypes$end<-as.integer(as.character(mendel_genotypes$end))
genotypes<-left_join(genotypes, mendel_genotypes,  by = c("chrom", "start", "end"))
return(genotypes)
}

##Histograms if/when needed
#trio_1<-data.frame(mendel$genotypes.chrom, mendel$genotypes.start, mendel$genotypes.end, mendel$genotypes.HG00512, mendel$genotypes.HG00513,
#                   mendel$genotypes.HG00514, mendel$HG00514_cons)
#trio_1$sample<-'HG00514'
#colnames(trio_1)<-c('chrom','start', 'end', 'father', 'mother', 'child', 'mendel_cons', 'sample')

#trio_2<-data.frame(mendel$genotypes.chrom, mendel$genotypes.start, mendel$genotypes.end, mendel$genotypes.HG00731, mendel$genotypes.HG00732,
#                   mendel$genotypes.HG00733, mendel$HG00733_cons)
#trio_2$sample<-'HG00733'
#colnames(trio_2)<-c('chrom','start', 'end', 'father', 'mother', 'child', 'mendel_cons', 'sample')

#trio_3<-data.frame(mendel$genotypes.chrom, mendel$genotypes.start, mendel$genotypes.end, mendel$genotypes.NA19238, mendel$genotypes.NA19239,
#                   mendel$genotypes.NA19240, mendel$NA19240_cons)
#trio_3$sample<-'NA19240'
#colnames(trio_3)<-c('chrom','start', 'end', 'father', 'mother', 'child', 'mendel_cons', 'sample')
#trio1_2<-rbind(trio_1,trio_2)
#all_trios<-rbind(trio1_2,trio_3)
##exclude 'noreads' regions
#all_trios_filter<-data.frame(all_trios[all_trios$mendel_cons!='noreads',])
##check whether the genotype is 'simple' or 'complex'
#all_trios_filter$call<-ifelse(((all_trios_filter$father=='1010' | all_trios_filter$father=='0101'| all_trios_filter$father=='1001'|  all_trios_filter$father=='0110') &
#                                 (all_trios_filter$mother=='1010' | all_trios_filter$mother=='0101'| all_trios_filter$mother=='1001'|  all_trios_filter$mother=='0110')&
#                                 (all_trios_filter$child=='1010' | all_trios_filter$child=='0101'| all_trios_filter$child=='1001'|  all_trios_filter$child=='0110')),
#                              'Simple', 'Complex')
#all_trios_filter$call_mendel<-NA
##check whether the inversion passed or failed the check
#all_trios_filter[all_trios_filter$call=='Simple',]$call_mendel<-ifelse(all_trios_filter[all_trios_filter$call=='Simple',]$mendel_cons=='YES',
#                                                                       'Simple_Pass', 'Simple_fail') 
#all_trios_filter[all_trios_filter$call=='Complex',]$call_mendel<-ifelse(all_trios_filter[all_trios_filter$call=='Complex',]$mendel_cons=='YES',
#                                                                        'Complex_Pass', 'Complex_fail') 

#library(ggplot2)
##theme_set(theme_classic())
#g <- ggplot(all_trios_filter, aes(call)) + scale_fill_brewer(palette = "Spectral")
#g + geom_histogram(aes(fill=call_mendel), 
#                   size=.1,
#                   binwidth = .1,
#                   stat = "count",
#                   col="black"
#)+   theme(axis.text.x = element_text(size =17),axis.text.y = element_text(size =17),  text = element_text(size=17)) +
#  labs(title="Mendelian_Test", x='Genotype', y='Count')+  theme(legend.title=element_blank())













