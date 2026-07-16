rm(list=ls())

library(HiCcompare)
library(dplyr)
library(data.table)

binsize <- "250kb"
comparisons <- list(c("alpha", "beta"), c("alpha", "STM"))

print(comparisons)

if (binsize == "1Mb") {
  A_min = 25
} else if (binsize == "500kb") {
  A_min = 12.5
} else if (binsize == "250kb") {
  A_min = 6.25
}

basedir <- "/home/anton/Bureau/PORE-C_repo/"
rdir <- paste0(basedir, "data/regions/")
cmdir <- paste0(basedir, "data/contact_maps/")


chromosomes <- c("NC_047559.1", 
                 "NC_047560.1", 
                 "NC_047561.1", 
                 "NC_047562.1", 
                 "NC_047563.1", 
                 "NC_047564.1", 
                 "NC_047565.1", 
                 "NC_047566.1", 
                 "NC_047567.1", 
                 "NC_047568.1")


for (comparison in comparisons) {
  
  sample1 <- comparison[1]
  sample2 <- comparison[2]

  print(paste("comparison", comparison))
  cat(paste(sample1, "vs", sample2))
  
  cm1 <- cooler2bedpe(path=paste0(cmdir, sample1, "/", sample1, "_", binsize, ".cool"))
  cm2 <- cooler2bedpe(path=paste0(cmdir, sample2, "/", sample2, "_", binsize, ".cool"))

  # create, filter and normalize comparison matrices
  for (chromosome in chromosomes) {
    
    sample1_chr <- cm1$cis[[chromosome]]
    sample2_chr <- cm2$cis[[chromosome]]
    
    chr_table <- create.hic.table(sample1_chr, sample2_chr, scale = TRUE)
    chr_table %>% filter(IF1 > 1) %>% filter(IF2 > 1)
    filter_params(chr_table)
    chr_table <- hic_loess(chr_table, Plot = T, Plot.smooth = T)
    table_chr <- hic_compare(chr_table, adjust.dist = TRUE, A.min = A_min, p.method = 'fdr', Plot = TRUE)

    table_name <- paste0("table_", chromosome)
    assign(table_name, table_chr)

  }
 
  tableList <- list()
  signifTableList <- list()
    
  for (chromosome in chromosomes) {
	  chromTable <- get(paste0("table_", chromosome))
	  tableList[[chromosome]] <- chromTable
	  signifTableList[[chromosome]] <- chromTable[chromTable$p.adj < 0.05, ]
  }
    
  all_bins <- do.call(rbind, tableList)
  all_signif <- do.call(rbind, signifTableList)
  
  compStr = paste(comparison[1], "vs", comparison[2])
  # coordinates corresponding to the pairs with significant interaction differences (to intersect with the genes)
  all_signif_pos1 <- all_signif %>% mutate(log_p = -log10(p.adj)) %>% select(chr1, start1, end1, A, log_p, M) %>% rename(chr = "chr1", start = "start1", end = "end1") %>% mutate(A = 0)
  all_signif_pos2 <- all_signif %>% mutate(log_p = -log10(p.adj)) %>% select(chr2, start2, end2, A, log_p, M) %>% rename(chr = "chr2", start = "start2", end = "end2") %>% mutate(A = 0)
  all_signif_bins <- rbind(all_signif_pos1, all_signif_pos2) %>% mutate()
    
  write.table(all_bins, file = file.path(rdir, paste0("all_bins_", sample1, "_", sample2, "_", binsize, ".bed")), 
              quote=F, row.names=F, col.names=T, sep="\t")
  write.table(all_signif_bins, file = file.path(rdir, paste0("all_signif_regions_", sample1, "_", sample2, "_", binsize, ".bed")), 
                quote=F, row.names=F, col.names=F, sep="\t")
    
  }
