rm(list=ls())

library(HiCcompare)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggpubr)

binsizes <- c("1Mb", "500kb", "250kb")
comparisons <- list(c("alpha", "beta"), c("alpha", "STM"))

print(comparisons)

basedir <- "/home/anton/Bureau/PORE-C_repo_Sanger/"
rdir <- paste0(basedir, "data/regions/")
cmdir <- paste0(basedir, "data/contact_maps/")
pdir <- paste0(basedir, "plots/HiCcompare")

chromosomes <- c("NC_088853.1",
		 "NC_088854.1",
		 "NC_088855.1",
		 "NC_088856.1",
		 "NC_088857.1",
		 "NC_088858.1",
		 "NC_088859.1",
		 "NC_088860.1",
		 "NC_088861.1",
		 "NC_088862.1")

for (binsize in binsizes) {

	for (comparison in comparisons) {

		if (binsize == "1Mb") {
		  A_min = 30
		} else if (binsize == "500kb") {
		  A_min = 15
		} else if (binsize == "250kb") {
		  A_min = 7.5
		}

	 
	  sample1 <- comparison[1]
	  sample2 <- comparison[2]
	  
	  pdir <- paste0(basedir, "plots/HiCcompare/", sample1, "_", sample2, "/", binsize, "_bins")
	  
	  print(paste("comparison : ", comparison))
	  cat(paste(sample1, "vs", sample2, "\n"))
	  
	  cm1 <- cooler2bedpe(path=paste0(cmdir, sample1, "/", sample1, "_", binsize, ".cool"))
	  cm2 <- cooler2bedpe(path=paste0(cmdir, sample2, "/", sample2, "_", binsize, ".cool"))

	  # create, filter and normalize comparison matrices
	  
	  pdf(file = file.path(pdir, paste0(sample1, "_", sample2, "_", binsize, "_report.pdf")))
	  
	  for (chromosome in chromosomes) {
	    
	    sample1_chr <- cm1$cis[[chromosome]]
	    sample2_chr <- cm2$cis[[chromosome]]
	    
	    chr_table <- create.hic.table(sample1_chr, sample2_chr, scale = TRUE)
	    chr_table %>% filter(IF1 > 1) %>% filter(IF2 > 1)
	    filter_params(chr_table)
	    chr_table <- hic_loess(chr_table, Plot = T, Plot.smooth = T)
	    # get the nth percentile of the A value
	    Amins = quantile(chr_table$A, probs = c(0.75))
	    print(Amins)
	    Amin = Amins[[1]]
	    table_chr <- hic_compare(chr_table, adjust.dist = TRUE, A.min = Amin, p.method = 'fdr', Plot = TRUE)

	    table_name <- paste0("table_", chromosome)
	    assign(table_name, table_chr)

	  }
	  
	  dev.off()
	  
	  tableList <- list()
	  signifTableList <- list()
	    
	  for (chromosome in chromosomes) {
		  chromTable <- get(paste0("table_", chromosome))
		  tableList[[chromosome]] <- chromTable
		  signifTableList[[chromosome]] <- chromTable[chromTable$p.adj < 0.05, ]
	  }
	    
	  all_bins <- do.call(rbind, tableList)
	  all_signif <- do.call(rbind, signifTableList)
	  
	  # distribution of A values
	  quant <- data.frame(quantile(all_bins$A, probs=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)))
	  colnames(quant) <-  c("A")
	  quant$pct = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
	  
	  pct <- ggplot(data = quant, aes(x=pct, y=A)) +
	    geom_col(fill="lightgrey",  linewidth = 0.5, fill="lightgrey", color="black")  +
	    theme_bw(base_size = 18) +
	    theme(panel.grid = element_blank(),
	          axis.text = element_text(color="black"),
	          axis.ticks = element_line(color="black"),
	          aspect.ratio = 1.2) +
	    scale_y_log10() +
	    labs(x = "percentile", y = "sum of contact frequencies (A value)")
	  
	  hist <- ggplot(data = all_bins, aes(y=A)) + 
	    geom_histogram(bins=30, linewidth = 0.5, fill="lightgrey", color="black") +
	    coord_flip() +
	    theme_bw(base_size = 18) +
	    theme(panel.grid = element_blank(),
	          axis.text = element_text(color="black"),
	          axis.ticks = element_line(color="black"),
	          aspect.ratio = 1.2) +
	    scale_x_log10() +
	    scale_y_log10() +
	    labs(x = "count", y="sum of contact frequencies (A value)")
	  
	  comb <- ggarrange(pct, hist)
	  ggsave(file.path(pdir, paste0("A_distribution_", sample1, "_", sample2, "_", binsize, ".pdf")),
	          plot = comb, device = cairo_pdf, width = 12, height = 8)
	  
	  
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
}
