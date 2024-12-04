library(magrittr)
library(GenomicRanges)

input <- snakemake@input[[1]]
output <- snakemake@output[[1]]
min <- snakemake@params[[1]]

gr <- readRDS(input)

gr2 <- gr %>%
  .[(gr %>%
      mcols %>%
      as.matrix %>%
      rowSums) >= min]

df <- data.frame(seqnames=seqnames(gr2),
                 starts=start(gr2)-1,
                 ends=end(gr2),
                 names=c(rep(".", length(gr2))),
                 scores=c(rep(".", length(gr2))),
                 strands=strand(gr2))

write.table(df, file=output, quote=F, sep="\t", row.names=F, col.names=F)
