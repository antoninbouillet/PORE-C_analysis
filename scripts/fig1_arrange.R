library(cowplot)
library(ggplot2)
library(magick)

rm(list = ls())

baseDir <- "/home/anton/Bureau/PORE-C_repo/plots/"
binsize <- "1Mb"
samples <- c("alpha", "beta", "STM")


readPdf <- function(pdfPath, page = 1, density = 600) {
	pdf <- image_read_pdf(pdfPath, page = page, density = density)
	pdf <- image_trim(pdf)
	return(as.raster(pdf))
}

# load contact vs distance plots
for (sample in samples) {

	cvd <- grid::rasterGrob(readPdf(paste0(baseDir, "contacts_distance/", sample, "/", sample, "_", binsize, "_cvd.pdf")), interpolate = F)
	assign(paste0(sample, "_cvd"), cvd)
}

# load an example of a contact map on chr3 : contact maps alpha + significant bins + difference alpha / beta and alpha / STM

boxplotChr3 <- grid::rasterGrob(readPdf(paste0(baseDir, paste0("contacts_distance/cvd_boxplot_", binsize, ".pdf")), page = 3), interpolate = F)
mapChr3 <- grid::rasterGrob(readPdf(paste0(baseDir, paste0("fanc_analysis/alpha/alpha_", binsize, "_insulation.pdf")), page = 3), interpolate = F)
mirrorChr3 <- grid::rasterGrob(readPdf(paste0(baseDir, paste0("fanc_comparison/mirror_zoom_chr3_", binsize, ".pdf"))), interpolate = F)

cvdY = 0.58
cvdW = 0.33
cvdScale = 0.9

fig1 <- ggdraw() +
  draw_plot(alpha_cvd, x = 0,   y = cvdY, width = cvdW, height = 0.5, scale = cvdScale) +
  draw_plot(beta_cvd, x = 0.33, y = cvdY, width = cvdW, height = 0.5, scale = cvdScale) +
  draw_plot(STM_cvd, x = 0.66,   y = cvdY, width = cvdW, height = 0.5, scale = cvdScale) +
  draw_plot(mapChr3, x = 0, y = 0.17, width = 0.5, height = 0.5, scale = 1) +
  draw_plot(boxplotChr3, x = 0.5, y = 0.3, width = 0.5, height = 0.5, scale = 0.8) +
  draw_plot(mirrorChr3, x = 0.5, y = -0.05, width = 0.5, height = 0.5, scale = 0.8) +

  draw_plot_label("A", x = 0.01, y = 0.98, size = 80, hjust = 0) +
  draw_plot_label("B", x = 0.01, y = 0.68, size = 80, hjust = 0) +
  draw_plot_label("C", x = 0.5, y = 0.68, size = 80, hjust = 0) +
  draw_plot_label("D", x = 0.5, y = 0.4, size = 80, hjust = 0)

ggsave2(file.path(baseDir, paste0("fig1_", binsize, ".pdf")),
         plot = fig1, device = cairo_pdf, width = 45, height = 45)
