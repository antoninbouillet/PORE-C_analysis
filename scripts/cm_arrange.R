library(cowplot)
library(ggplot2)
library(magick)

rm(list = ls())

baseDir <- "/home/anton/Bureau/PORE-C_repo_Sanger/plots/contact_maps"
binsizes <- c("1Mb", "500kb", "250kb")
samples <- c("alpha", "beta", "STM")


readPdf <- function(pdfPath, page = 1, density = 250) {
	pdf <- image_read_pdf(pdfPath, page = page, density = density)
	pdf <- image_trim(pdf)
	return(as.raster(pdf))
}

for (binsize in binsizes) {

	for (sample in samples) {
		
		for (i in seq(10)) {
			cm <- grid::rasterGrob(readPdf(paste0(baseDir, "/", sample, "_",  binsize, ".pdf"), page = i), interpolate = F)
			assign(paste0(sample, "_chr", i), cm)
		}

	}
	x0 = 0
	x1 = 0.12
	x2 = 0.24

	allCm <- ggdraw() +

		draw_plot(alpha_chr1, x = x0, y = 0.9, width = 0.33, height = 0.1) +
		draw_plot(beta_chr1, x = x1, y = 0.9, width = 0.33, height = 0.1) +
		draw_plot(STM_chr1, x = x2, y = 0.9, width = 0.33, height = 0.1) +

		draw_plot(alpha_chr2, x = x0, y = 0.8, width = 0.33, height = 0.1) +
		draw_plot(beta_chr2, x = x1, y = 0.8, width = 0.33, height = 0.1) +
		draw_plot(STM_chr2, x = x2, y = 0.8, width = 0.33, height = 0.1) +

		draw_plot(alpha_chr3, x = x0, y = 0.7, width = 0.33, height = 0.1) +
		draw_plot(beta_chr3, x = x1, y = 0.7, width = 0.33, height = 0.1) +
		draw_plot(STM_chr3, x = x2, y = 0.7, width = 0.33, height = 0.1) +

		draw_plot(alpha_chr4, x = x0, y = 0.6, width = 0.33, height = 0.1) +
		draw_plot(beta_chr4, x = x1, y = 0.6, width = 0.33, height = 0.1) +
		draw_plot(STM_chr4, x = x2, y = 0.6, width = 0.33, height = 0.1) +

		draw_plot(alpha_chr5, x = x0, y = 0.5, width = 0.33, height = 0.1) +
		draw_plot(beta_chr5, x = x1, y = 0.5, width = 0.33, height = 0.1) +
		draw_plot(STM_chr5, x = x2, y = 0.5, width = 0.33, height = 0.1) +

		draw_plot(alpha_chr6, x = x0, y = 0.4, width = 0.33, height = 0.1) +
		draw_plot(beta_chr6, x = x1, y = 0.4, width = 0.33, height = 0.1) +
		draw_plot(STM_chr6, x = x2, y = 0.4, width = 0.33, height = 0.1) +

		draw_plot(alpha_chr7, x = x0, y = 0.3, width = 0.33, height = 0.1) +
		draw_plot(beta_chr7, x = x1, y = 0.3, width = 0.33, height = 0.1) +
		draw_plot(STM_chr7, x = x2, y = 0.3, width = 0.33, height = 0.1) +

		draw_plot(alpha_chr8, x = x0, y = 0.2, width = 0.33, height = 0.1) +
		draw_plot(beta_chr8, x = x1, y = 0.2, width = 0.33, height = 0.1) +
		draw_plot(STM_chr8, x = x2, y = 0.2, width = 0.33, height = 0.1) +

		draw_plot(alpha_chr9, x = x0, y = 0.1, width = 0.33, height = 0.1) +
		draw_plot(beta_chr9, x = x1, y = 0.1, width = 0.33, height = 0.1) +
		draw_plot(STM_chr9, x = x2, y = 0.1, width = 0.33, height = 0.1) +

		draw_plot(alpha_chr10, x = x0, y = 0, width = 0.33, height = 0.1) +
		draw_plot(beta_chr10, x = x1, y = 0, width = 0.33, height = 0.1) +
		draw_plot(STM_chr10, x = x2, y = 0, width = 0.33, height = 0.1)

	ggsave2(file.path(baseDir, paste0("contact_maps", binsize, ".pdf")),
		 plot = allCm, device = cairo_pdf, width = 45, height = 45)


}


