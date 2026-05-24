rm(list = ls())

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggnewscale)
library(ggvenn)

binsizes <- c("1Mb", "500kb", "250kb")
comparisons <- list(c("alpha", "beta"), c("alpha", "STM"))
color1 <- "#5948d9"
colorBeta <- "#48D959"
colorSTM <- "#D95948"

zoomChr = "NC_047561.1"
zoomStart = 2e7
zoomEnd = 5.5e7

baseDir <- "/home/anton/Bureau/PORE-C_repo/"
rDir <- paste0(baseDir, "data/regions")

for (binsize in binsizes) {

  if (binsize == "1Mb") {
    pixelsize <- 0.47
    zoomPixelsize <- 4.7
  } else if (binsize == "500kb") {
    pixelsize <- 0.09
    zoomPixelsize <- 2
  } else if (binsize == "250kb") {
    pixelsize <- 0.01
    zoomPixelsize <- 0.75
  }
  
  for (comparison in comparisons) {
  
  	sample1 <- comparison[1]
  	sample2 <- comparison[2]
  	
  	pdir <- paste0(baseDir, "/plots/HiCcompare/", sample1, "_", sample2, "/", binsize, "_bins")
  	
  	status1 <- paste0(sample1, "+")
  	status2 <- paste0(sample2, "+")
  
  	label1 <- "α"
  	if (sample2 == "beta") {
  		label2 <- "β"
  		color2 <- colorBeta
  		levels <- c("α+", "β+")
  	} else {
  		label2 <- "STM"
  		color2 <- colorSTM
  		levels <- c("α+", "STM+")
  	}
  
  	all_bins <- fread(file=file.path(rDir, paste0("all_bins_", sample1, "_", sample2, "_", binsize, ".bed")), header=T) %>% 
  	  mutate(status = ifelse(M<0, status1, status2),
  		 is_signif = ifelse(p.adj < 0.05, TRUE, FALSE),
  		 pval_log = -log10(p.adj))
  
  	all_bins$status <- sub("alpha\\+", "α+",  all_bins$status)
  	if (sample2 == "beta") {
  		all_bins$status <- sub("beta\\+", "β+", all_bins$status)
  	}
  
  
  	all_bins$status <- factor(all_bins$status, levels = levels)
  	
  	# save bins with signifiant interaction differences to check overlap between comparisons
  	signif <- all_bins[all_bins$is_signif == T, ]
  	signifName <- paste0("signif_", sample1, "_", sample2)
	  assign(signifName, signif)

  	volcanoplot <- ggplot() + 
  	  geom_point(data=signif,
  		     aes(x = adj.M, y = pval_log, color = status),
  		     size=1, alpha=0.7) +
  	  geom_point(data=all_bins[all_bins$pval_log > 0 & all_bins$is_signif == F, ],
  		     aes(x = adj.M, y = pval_log),
  		     size=1, alpha=0.7, color="grey") +
  	  geom_hline(yintercept = -log10(0.05), linetype="dashed", linewidth=0.2) +
  	  theme_bw(base_size=22) +
  	  theme(aspect.ratio = 1,
  		panel.grid = element_line(linetype="dashed", linewidth=0.2),
  		legend.title = element_blank()) +
  	  scale_y_continuous(limits=c(0, 1.1*max(all_bins$pval_log)), expand=c(0,0)) +
  	  scale_x_continuous(limits=c(-3,3)) +
  	  scale_color_manual(values=c(color1, color2)) +
  	  labs(y = "-log10(P-adj)", x = "adjuted M value")
  	
  	  ggsave(file.path(pdir, paste0("volcanoplot_", sample1, "_", sample2, "_", binsize, ".pdf")), 
  	       plot = volcanoplot, device = cairo_pdf, width = 8, height = 6)

	  # also save the plots as RDS objects to arrange them later
	  saveRDS(volcanoplot, file = file.path(pdir, paste0("volcanoplot_", sample1, "_", sample2, "_", binsize, ".rds")))
  
  	  maps <- ggplot(data = all_bins)  + 
  	  theme_bw(base_size=7) +
  	  theme(aspect.ratio = 1, panel.grid=element_blank(),
  		strip.background = element_rect(fill="white", color="white"),
  		axis.text = element_text(color="black"),
  		panel.border = element_rect(color="black"),
  		axis.ticks = element_line(color="black")) +
      geom_point(aes(x=NA, y=NA, color=status)) +
      scale_color_manual(name="", values=c(color1, color2)) +
  	  guides(color=guide_legend(order=1)) +
      new_scale_color() +
  	  geom_point(aes(x=start1, y=start2, color=adj.M), shape = 15, size=pixelsize) + 
  	  scale_color_gradient2(low=color1, mid="white", high=color2) +
      new_scale_color() +
  	  geom_point(aes(x=start2, y=start1, color=pval_log), shape = 15, size=pixelsize) +
  	  scale_color_gradient(name="-log10(P-adj)", low="white", high="black") +
  	  labs(x="", y="") +
  	  facet_wrap(~chr1, scale = "free") +
  	  scale_y_continuous(expand=c(-1.015,1.015)) +
  	  scale_x_continuous(expand=c(-1.015,1.015))
  
  	ggsave(file.path(pdir, paste0("maps_HiCcompare_", sample1, "_", sample2, "_", binsize, ".pdf")), 
  	       plot = maps, device = cairo_pdf, width = 8, height = 6)
  	
  	# zoom on one region of chr3 as an example
  	mapZoom <- ggplot(data = all_bins[all_bins$chr1 == zoomChr & 
  	                                 all_bins$start1 > zoomStart & 
  	                                 all_bins$start1 < zoomEnd &
  	                                 all_bins$start2 > zoomStart & 
  	                                 all_bins$start2 < zoomEnd, ])  + 
  	  theme_bw(base_size=18) +
  	  theme(aspect.ratio = 1, panel.grid=element_blank(),
  	        strip.background = element_rect(fill="white", color="white"),
  	        axis.text = element_text(color="black"),
  	        panel.border = element_rect(color="black"),
  	        axis.ticks = element_line(color="black")) +
  	  geom_point(aes(x=NA, y=NA, color=status)) +
  	  scale_color_manual(name="", values=c(color1, color2)) +
  	  guides(color=guide_legend(order=1)) +
  	  new_scale_color() +
  	  geom_point(aes(x=start1, y=start2, color=adj.M), shape = 15, size=zoomPixelsize) + 
  	  scale_color_gradient2(low=color1, mid="white", high=color2) +
  	  new_scale_color() +
  	  geom_point(aes(x=start2, y=start1, color=pval_log), shape = 15, size=zoomPixelsize) +
  	  scale_color_gradient(name="-log10(P-adj)", low="white", high="black") +
  	  labs(x="", y="") +
  	  scale_y_continuous(expand=c(-1.015,1.015)) +
  	  scale_x_continuous(expand=c(-1.015,1.015))
  	
  	ggsave(file.path(pdir, paste0("map_zoom_HiCcompare_", sample1, "_", sample2, "_", binsize, ".pdf")), 
  	       plot = mapZoom, device = cairo_pdf, width = 8, height = 6)

	saveRDS(mapZoom, file = file.path(pdir, paste0("map_zoom_HiCcompare_", sample1, "_", sample2, "_", binsize, ".rds")))
  	
  	mean_IF1 <- all_bins %>% group_by(chr1, D) %>% select(chr1, D, adj.IF1) %>% summarize(meanIF=mean_se(adj.IF1)) %>% mutate(ech=label1)
  	mean_IF2 <- all_bins %>% group_by(chr1, D) %>% select(chr1, D, adj.IF2) %>% summarize(meanIF=mean_se(adj.IF2)) %>% mutate(ech=label2)
  	all_IF <- rbind(mean_IF1, mean_IF2)
  	all_IF$ech <- factor(all_IF$ech, levels=c(label1, label2))
  
  	# bin the distance and compare the means of interaction frequency for each bin
  	all_IF <- all_IF %>% mutate(bins = ntile(D, 5))
  	all_IF$bins <- as.factor(all_IF$bins)
  	
  	cdbin <- ggplot(data = all_IF, aes(x = bins, y = meanIF$y, fill=ech)) +
  	  geom_boxplot(outliers = F, linewidth=0.25) +
  	  theme_bw(base_size=10) +
  	  # stat_compare_means(size=1, angle=90) +
  	  theme(strip.background = element_rect(fill="white", color="white"),
  	        axis.text = element_text(color="black"),
  	        panel.border = element_rect(color="black"),
  	        axis.ticks = element_line(color="black"),
  	        aspect.ratio = 1,
  	        panel.grid = element_blank(),
  	        legend.title = element_blank()) +
  	  labs(x="distance (binned)", y="mean interaction frequency") +
  	  facet_wrap(~chr1, scale = "free_y") +
  	  scale_y_log10() +
  	  scale_fill_manual(values=c(color1, color2))
  	
  	ggsave(file.path(pdir, paste0("contact_distance_binned_HiCcompare_", sample1, "_", sample2, "_", binsize, ".pdf")), 
  	       plot = cdbin, device = cairo_pdf, width = 8, height = 6)

  	cd <- ggplot(data=all_IF) +
  	  theme_bw() +
  	  theme(strip.background = element_rect(fill="white", color="white"),
  		axis.text = element_text(color="black"),
  		panel.border = element_rect(color="black"),
  		axis.ticks = element_line(color="black"),
  		aspect.ratio = 1,
  		panel.grid = element_blank(),
  		legend.title = element_blank()) +
  	  geom_line(aes(x=D+1, y=meanIF$y, color=ech), linewidth=0.1) +
  	  geom_line(aes(x=D+1, y=meanIF$ymin, color=ech), linetype="dashed", linewidth=0.1) +
  	  geom_line(aes(x=D+1, y=meanIF$ymax, color=ech),linetype="dashed", linewidth=0.1) +
  	  facet_wrap(~chr1, scale="free_x") +
  	  labs(x=paste0("Distance (", binsize, " bins)"), y="Interaction frequency") +
  	  scale_y_log10() +
  	  scale_x_log10() +
  	  scale_color_manual(values=c(color1, color2))
  
  	ggsave(file.path(pdir, paste0("contact_distance_HiCcompare_", sample1, "_", sample2, "_", binsize, ".pdf")), 
  	       plot = cd, device = cairo_pdf, width = 8, height = 6)

  	# filtering the HIC tables by A value (mean of IF1 and IF2) tends to favor pairs with
  	# low IF (which are also the most distant) in the matrix with the lowest overall read counts (beta / STM)
  
  	comp_d <- ggplot(data=all_bins[all_bins$p.adj < 0.05, ], aes(x=status, y=D, fill=status)) +
  	  geom_boxplot(width=0.4, color="black", linewidth=0.8, outliers=F) +
  	  theme_bw(base_size=22) +
  	  theme(panel.grid=element_blank(),
  		panel.border=element_rect(color="black"),
  		axis.text = element_text(color="black"),
  		axis.ticks = element_line(color="black"),
  		aspect.ratio=1.2,
  		legend.position="None") +
  	  labs(y=paste0("Distance (", binsize, " bins)"), x="") +
  	  stat_compare_means(size=6) +
  	  scale_fill_manual(values=c(color1, color2))
  	
  	ggsave(file.path(pdir, paste0("comp_distance_", sample1, "_", sample2, "_", binsize, ".pdf")), 
  	       plot = comp_d, device = cairo_pdf, width = 8, height = 6)

	saveRDS(comp_d, file = file.path(pdir, paste0("comp_distance_", sample1, "_", sample2, "_", binsize, ".rds")))
  }
  
  # venn diagrams
  
  pdir <- paste0(baseDir, "plots/HiCcompare/overlap")
  comp <- full_join(signif_alpha_beta, signif_alpha_STM,
                    by=join_by(chr1, start1, start2),
                    suffix = c(".alpha_beta", ".alpha_STM")) 
  
  comp <- comp %>% 
    mutate(

      # total overlap
      alpha_beta = !is.na(status.alpha_beta),
      alpha_STM = !is.na(status.alpha_STM),

      # alpha +  
      alpha_beta_alpha = !is.na(status.alpha_beta) & status.alpha_beta == "α+",
      alpha_STM_alpha  = !is.na(status.alpha_STM) & status.alpha_STM == "α+",
      
      # beta + / STM + 
      alpha_beta_beta = !is.na(status.alpha_beta) & status.alpha_beta == "β+",
      alpha_STM_STM   = !is.na(status.alpha_STM) & status.alpha_STM == "STM+",

    )
  
    nb_size = 6
    label_size = 8

    p_all <- ggplot(data=comp) +
    geom_venn(aes(A = alpha_beta,
                  B = alpha_STM),
        text_size = nb_size,
	      set_name_size = 0) +
    coord_fixed() +
    labs(title="") +
    theme_void(base_size=10) +
    annotate("text", x=-0.8, y=1.3, label="α vs β", size=label_size) +
    annotate("text", x=0.8, y=1.3, label="α vs STM", size=label_size)

  
  ggsave(file.path(pdir, paste0("Venn_total", binsize, ".pdf")),
         plot = p_all, device = cairo_pdf, width = 8, height = 6)

  saveRDS(p_all, file = file.path(pdir, paste0("Venn_total", binsize, ".rds")))

  p_alpha <- ggplot(data=comp[comp$status.alpha_beta == "α+" |
                                comp$status.alpha_STM == "α+", ]) + 
    geom_venn(aes(A = alpha_beta_alpha,
                  B = alpha_STM_alpha),
	      fill_color = c(color1, color1),
	      text_size = nb_size,
	      set_name_size = 0)+
    coord_fixed() +
    theme_void(base_size=10) +
    annotate("text", x=-0.9, y=1.3, label="α vs β (α+)", size=label_size) +
    annotate("text", x=0.9, y=1.3, label="α vs STM (α+)", size=label_size)
  
  ggsave(file.path(pdir, paste0("Venn_alpha_alpha_", binsize, ".pdf")),
         plot = p_alpha, device = cairo_pdf, width = 8, height = 6)

  saveRDS(p_alpha, file = file.path(pdir, paste0("Venn_alpha_alpha_", binsize, ".rds")))
  
  p_other <- ggplot(data=comp[comp$status.alpha_beta == "β+" |
                                comp$status.alpha_STM == "STM+", ]) + 
    geom_venn(aes(A = alpha_beta_beta, 
                  B = alpha_STM_STM),
              fill_color = c(colorBeta, colorSTM),
	      text_size = nb_size,
	      set_name_size = 0) +
    coord_fixed() +
    theme_void(base_size=10) +
    theme(plot.title = element_text(size = 10)) +
    annotate("text", x=-0.9, y=1.3, label="α vs β (β+)", size=label_size) +
    annotate("text", x=0.9, y=1.3, label="α vs STM (STM+)", size=label_size)
  
  ggsave(file.path(pdir, paste0("Venn_beta_STM_", binsize, ".pdf")),
         plot = p_other, device = cairo_pdf, width = 8, height = 6)

  saveRDS(p_other, file = file.path(pdir, paste0("Venn_beta_STM_", binsize, ".rds")))
}
