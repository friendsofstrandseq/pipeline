log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

library(data.table)
library(ggplot2)
library(gridExtra)

br.probs <- Reduce(rbind, lapply(snakemake@input[["breakpoint_ll"]], function(x) fread(paste("zcat", x))))

br.ll.hist <- ggplot(br.probs, aes(x=log(br_ll), fill=breakpoint_type!="no_br"))+geom_histogram(aes(y=..density..))
br.ll.density <- ggplot(br.probs, aes(x=log(br_ll), col=breakpoint_type!="no_br"))+geom_density()

ggsave(snakemake@output[["plot_breakpoint_ll"]], grid.arrange(br.ll.hist, br.ll.density, nrow=2,  ncol=1))

