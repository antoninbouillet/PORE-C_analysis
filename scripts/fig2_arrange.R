library(cowplot)
library(egg)
library(ggplot2)
library(grid)

rm(list=ls())

# Figure 2 : volcanoplot -> contact map -> GO -> Venn diagram for alpha vs beta / alpha vs STM

# load plots

baseDir <- "/home/anton/Bureau/PORE-C_repo/plots/"
binsizes <- c("1Mb", "500kb", "250kb")
pDir <- paste0(baseDir, "HiCcompare/")


chrSize = 1

zoomChr = "NC_047561.1"
zoomStart = 2e7
zoomEnd = 5.5e7

chrDiagram <- ggplot() +
  ## Chromosome
  geom_segment(
    aes(x = 0, y = 0,
        xend = 58.319e6, yend = 0),,
    linewidth = 10,
    colour = "grey85",
    lineend = "round"
  ) +

  ## Highlighted region
  geom_segment(
    aes(x = zoomStart, y = 0,
        xend = zoomEnd, yend = 0),
    linewidth = 10,
    colour = "#4F81BD",
    lineend = "butt"
  ) +
  labs(x = "", y = "", title=zoomChr) +
  theme_classic(base_size = 24) +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    aspect.ratio = 0.065
  )

for (binsize in binsizes) {

	ABdir <- paste0(pDir, paste0("alpha_beta/", binsize, "_bins"))
	ASTMdir <- paste0(pDir, "alpha_STM/", binsize, "_bins")

	vSize = 1
	volcanoplotAB <- readRDS(file = file.path(ABdir, paste0("volcanoplot_alpha_beta_", binsize, ".rds")))
	volcanoplotASTM <- readRDS(file = file.path(ASTMdir, paste0("volcanoplot_alpha_STM_", binsize, ".rds")))

	dSize = 0.5
	distAB <- readRDS(file = file.path(ABdir, paste0("comp_distance_alpha_beta_", binsize, ".rds")))
	distASTM <- readRDS(file = file.path(ASTMdir, paste0("comp_distance_alpha_STM_", binsize, ".rds")))

	cSize = 1
	contactMapAB <- readRDS(file = file.path(ABdir, paste0("map_zoom_HiCcompare_alpha_beta_", binsize, ".rds")))
	contactMapASTM <- readRDS(file = file.path(ASTMdir, paste0("map_zoom_HiCcompare_alpha_STM_", binsize, ".rds")))

	gSize = 1
	goAB <- readRDS(file = file.path(ABdir, paste0("GO_terms_bins_", binsize, "_alpha_beta_all.rds")))
	goASTM <- readRDS(file = file.path(ASTMdir, paste0("GO_terms_bins_", binsize, "_alpha_STM_all.rds")))

	# build a representation of ch3 with karyoploteR


	#vennSize = 0.8
	#vennAA <- readRDS(file = file.path(paste0(pDir, "overlap"), paste0("Venn_alpha_alpha_", binsize, ".rds")))
	#vennBSTM <- readRDS(file = file.path(paste0(pDir, "overlap"), paste0("Venn_beta_STM_", binsize, ".rds")))

	sizeMod = 1.2
	fig2 <- ggdraw() +
	  
	  draw_plot(volcanoplotAB,    x = 0,   y = 0.7, width = 0.6 * sizeMod, height = 0.24 * sizeMod) +
	  draw_plot(volcanoplotASTM,  x = 0.3, y = 0.7, width = 0.5 * sizeMod, height = 0.24 * sizeMod) +
	  
	  draw_plot(distAB,           x = 0.225, y = 0.835, width = 0.22 * sizeMod, height = 0.11 * sizeMod) +
	  draw_plot(distASTM,         x = 0.455, y = 0.835, width = 0.22 * sizeMod, height = 0.11 * sizeMod) +
	  
	  draw_plot(chrDiagram,       x = 0.31, y = 0.57, width = 0.25 * sizeMod, height = 0.12 * sizeMod) +

	  draw_plot(contactMapAB,     x = 0.07, y = 0.315, width = 0.5 * sizeMod, height = 0.24 * sizeMod) +
	  draw_plot(contactMapASTM,   x = 0.3, y = 0.315, width = 0.5 * sizeMod, height = 0.24 * sizeMod) +
	  
	  draw_plot(goAB,             x = 0.06, y = 0, width = 0.5 * sizeMod, height = 0.24 * sizeMod) +
	  draw_plot(goASTM,           x = 0.292 , y = 0, width = 0.5 * sizeMod, height = 0.24 * sizeMod) +

	  # draw_plot(vennAA,          x = 0.094, y = 0.05, width = 0.33, height = 0.17) +
	  # draw_plot(vennBSTM,        x = 0.333, y = 0.05, width = 0.3, height = 0.17) +
	  
	  draw_plot_label("A", x = 0.22, y = 0.99, size = 48, hjust = 0) +
	  draw_plot_label("B", x = 0.22, y = 0.7, size = 48, hjust = 0) +
	  draw_plot_label("C", x = 0.22, y = 0.32, size = 48, hjust = 0)
	  # draw_plot_label("D", x = 0.13, y = 0.24, size = 48, hjust = 0)

	ggsave2(file.path(baseDir, paste0("fig2_", binsize, ".pdf")),
		 plot = fig2, device = cairo_pdf, width = 40, height = 25)
}
