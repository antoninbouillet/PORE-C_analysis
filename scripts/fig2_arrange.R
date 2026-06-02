library(cowplot)
library(egg)
library(ggplot2)

rm(list=ls())

# Figure 2 : volcanoplot -> contact map -> GO -> Venn diagram for alpha vs beta / alpha vs STM
# at a 500kb resolution

# load plots

baseDir <- "/home/anton/Bureau/PORE-C_repo/plots/"
pDir <- paste0(baseDir, "HiCcompare/")

ABdir <- paste0(pDir, "alpha_beta/500kb_bins")
ASTMdir <- paste0(pDir, "alpha_STM/500kb_bins")

vSize = 1
volcanoplotAB <- readRDS(file = file.path(ABdir, "volcanoplot_alpha_beta_500kb.rds"))
volcanoplotASTM <- readRDS(file = file.path(ASTMdir, "volcanoplot_alpha_STM_500kb.rds"))

dSize = 0.5
distAB <- readRDS(file = file.path(ABdir, "comp_distance_alpha_beta_500kb.rds"))
distASTM <- readRDS(file = file.path(ASTMdir, "comp_distance_alpha_STM_500kb.rds"))

cSize = 1
contactMapAB <- readRDS(file = file.path(ABdir, "map_zoom_HiCcompare_alpha_beta_500kb.rds"))
contactMapASTM <- readRDS(file = file.path(ASTMdir, "map_zoom_HiCcompare_alpha_STM_500kb.rds"))

gSize = 1
goAB <- readRDS(file = file.path(ABdir, "GO_terms_bins_500kb_alpha_beta_all.rds"))
goASTM <- readRDS(file = file.path(ASTMdir, "GO_terms_bins_500kb_alpha_STM_all.rds"))

vennSize = 0.8
vennAA <- readRDS(file = file.path(paste0(pDir, "overlap"), "Venn_alpha_alpha_500kb.rds"))
vennBSTM <- readRDS(file = file.path(paste0(pDir, "overlap"), "Venn_beta_STM_500kb.rds"))

scales <- c(vSize, vSize,
            cSize, cSize,
            gSize, gSize,
            vennSize, vennSize)

axes <- c("l", "l",
          "l", "l", 
          "l", "l"
          )

# using cowplot

#fig2 <- plot_grid(volcanoplotAB, volcanoplotASTM,
#		  contactMapAB, contactMapASTM, 
#		  goAB, goASTM,
#		  vennAA, vennBSTM,
#          labels = c("A", "", "B", "", "C", "", "D", ""),
#          ncol = 2, nrow = 4,
#		  # axis = axes,
#		  scale = scales,
#		  label_size = 48,
#		  hgap = 0
#		  )

fig2 <- ggdraw() +
  
  draw_plot(volcanoplotAB,    x = 0.006, y = 0.75, width = 0.5, height = 0.24) +
  draw_plot(volcanoplotASTM, x = 0.239, y = 0.75, width = 0.5, height = 0.24) +
  
  draw_plot(distAB, x = 0.135, y = 0.865, width = 0.22, height = 0.11) +
  draw_plot(distASTM, x = 0.360, y = 0.865, width = 0.22, height = 0.11) +
  
  draw_plot(contactMapAB,    x = 0.015, y = 0.50, width = 0.5, height = 0.24) +
  draw_plot(contactMapASTM, x = 0.238, y = 0.50, width = 0.5, height = 0.24) +
  
  draw_plot(goAB,            x = 0.002, y = 0.25, width = 0.5, height = 0.24) +
  draw_plot(goASTM,          x = 0.225, y = 0.25, width = 0.5, height = 0.24) +
  
  draw_plot(vennAA,          x = 0.094, y = 0.05, width = 0.33, height = 0.17) +
  draw_plot(vennBSTM,        x = 0.333, y = 0.05, width = 0.3, height = 0.17) +
  
  draw_plot_label("A", x = 0.13, y = 0.99, size = 48, hjust = 0) +
  draw_plot_label("B", x = 0.13, y = 0.74, size = 48, hjust = 0) +
  draw_plot_label("C", x = 0.13, y = 0.49, size = 48, hjust = 0) +
  draw_plot_label("D", x = 0.13, y = 0.24, size = 48, hjust = 0)

ggsave2(file.path(baseDir, "fig2.pdf"),
         plot = fig2, device = cairo_pdf, width = 40, height = 25)
