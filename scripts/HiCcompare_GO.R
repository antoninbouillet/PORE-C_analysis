#GO terms sur gènes situés dans régions à interaction différentielle

rm(list = ls())


library(ontologyIndex)
library(clusterProfiler)
library(data.table)
library(Cairo)
library(ggplot2)
library(dplyr)

binsize = "250kb"
baseDir <- "/home/anton/Bureau/PORE-C_repo/"
rDir <- paste0(baseDir, "data/regions")

comparisons <- list(c("alpha", "beta"), c("alpha", "STM"))

# load ontology data
ontology <- get_ontology("/home/anton/Bureau/MG_genome/go.obo")
gaf_data <- read.table("/home/anton/Bureau/MG_genome/GCF_902806645.1/GO/GCF_902806645.1_cgigas_uk_roslin_v1_gene_ontology.gaf", 
		       comment = "!", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

colnames(gaf_data) <- c("DB", "GeneID", "Symbol", "Qualifier", "GO_ID", "Reference", 
                         "Evidence_Code", "With_From", "Aspect", "Gene_Name", 
                         "Gene_Synonym", "Type", "Taxon", "Date", "Assigned_By", 
                         "Annot_Ext", "Gene_Product_Form_ID")
term2name <- data.frame(ontology$id, ontology$name) %>% rename(GO_ID = "ontology.id")
gaf_data <- gaf_data %>% filter(GO_ID %in% ontology$id)
gene_go_mapping <- gaf_data %>% left_join(term2name, by="GO_ID")
universe <- as.character(unique(gene_go_mapping$GeneID))


for (comparison in comparisons) {

	sample1 <- comparison[1]
	sample2 <- comparison[2]

	pdir <- paste0(baseDir, "plots/HiCcompare/", sample1, "_", sample2, "/", binsize, "_bins")

	status1 <- paste0(sample1, "+")
	status2 <- paste0(sample2, "+")

	label1 <- "α"
	if (sample2 == "beta") {
		label2 <- "β"
	} else {
		label2 <- sample2
	}

	all_signif_genes <- fread(file=file.path(rDir, paste0("all_signif_genes_", sample1, "_", sample2, "_", binsize, ".bed"))) %>% 
	  mutate(status = ifelse(V6<0, status1, status2))

	all_signif_1 <- all_signif_genes[all_signif_genes$status==status1, ]

	all_signif_2 <- all_signif_genes[all_signif_genes$status==status2, ]

	liste_genes_map <- gene_go_mapping %>% filter(Symbol %in% all_signif_genes$V11) 
	n_mapped <-  n_distinct(liste_genes_map$Symbol)
	n_total <- n_distinct(all_signif_genes$V11)

	liste_genes_enrich <-enricher(
	  as.character(liste_genes_map$GeneID),
	  pvalueCutoff = 0.1,
	  pAdjustMethod = "BH",
	  universe = universe,
	  minGSSize = 20,
	  maxGSSize = 38000,
	  qvalueCutoff = 0.1,
	  gson = NULL,
	  TERM2GENE = gene_go_mapping %>% select(GO_ID, GeneID),
	  TERM2NAME = gene_go_mapping %>% select(GO_ID, ontology.name))

	plotAll <- dotplot(liste_genes_enrich) +
	  theme_bw(base_size=14, base_rect_size=1) +
	  theme(panel.grid=element_blank(), 
	        panel.border = element_rect(color="black"),
	        axis.ticks = element_line(color="black"),
	        axis.text = element_text(color="black"),
		plot.title = element_text(size=12, hjust=0.5), 
		aspect.ratio  = 1.2) +
	  labs(title = paste(label1, "vs", label2, "(", binsize, "bins) :", n_mapped," / ", n_total, "genes assigned"))

	ggsave(file.path(pdir, paste0("GO_terms_bins_", binsize, "_", sample1, "_", sample2, "_all.pdf")), 
	       plot = plotAll, device = cairo_pdf, width = 8, height = 6)

	#genes sample1+ (alpha+)

	liste_genes_map_1 <- gene_go_mapping %>% filter(Symbol %in% all_signif_1$V11) 
	n_mapped_1 <- n_distinct(liste_genes_map_1$Symbol)
	n_total_1 <- n_distinct(all_signif_1$V11)


	liste_genes_enrich_1 <-enricher(
	  as.character(liste_genes_map_1$GeneID),
	  pvalueCutoff = 0.1,
	  pAdjustMethod = "BH",
	  universe = universe,
	  minGSSize = 20,
	  maxGSSize = 38000,
	  qvalueCutoff = 0.1,
	  gson = NULL,
	  TERM2GENE = gene_go_mapping %>% select(GO_ID, GeneID),
	  TERM2NAME = gene_go_mapping %>% select(GO_ID, ontology.name))

	plotS1 <- dotplot(liste_genes_enrich_1) +
	  theme_bw(base_size=14, base_rect_size=1) +
	  theme(panel.grid=element_blank(),
	        panel.border = element_rect(color="black"),
	        axis.ticks = element_line(color="black"),
	        axis.text = element_text(color="black"),
		plot.title = element_text(size=12, hjust=0.5), 
		aspect.ratio  = 1.2) +
	  labs(title = paste(label1, "vs", label2,"(", binsize, "bins", label1, "+ ) :", n_mapped," / ", n_total_1, "genes assigned"))

	ggsave(file.path(pdir, paste0("GO_terms_bins_", binsize, "_", sample1, "_", sample2, "_", status1, ".pdf")), 
	       plot = plotS1, device = cairo_pdf, width = 8, height = 6)


	#genes sample2+ (beta / STM)

	liste_genes_map_2 <- gene_go_mapping %>% filter(Symbol %in% all_signif_2$V11) 
	n_mapped_2 <- n_distinct(liste_genes_map_2$Symbol)
	n_total_2 <- n_distinct(all_signif_2$V11)

	liste_genes_enrich_2 <-enricher(
	  as.character(liste_genes_map_2$GeneID),
	  pvalueCutoff = 0.1,
	  pAdjustMethod = "BH",
	  universe = universe,
	  minGSSize = 20,
	  maxGSSize = 38000,
	  qvalueCutoff = 0.1,
	  gson = NULL,
	  TERM2GENE = gene_go_mapping %>% select(GO_ID, GeneID),
	  TERM2NAME = gene_go_mapping %>% select(GO_ID, ontology.name))

	plotS2 <- dotplot(liste_genes_enrich_2) +
	  theme_bw(base_size=14, base_rect_size=1) +
	  theme(panel.grid=element_blank(),
	        panel.border = element_rect(color="black"),
	        axis.ticks = element_line(color="black"),
	        axis.text = element_text(color="black"), 
		plot.title = element_text(size=12, hjust=0.5), 
		aspect.ratio  = 1.2) +
	  labs(title = paste(label1, "vs", label2, "(", binsize, "bins", label2, "+ ) :", n_mapped_2," / ", n_total_2, "genes assigned"))

	ggsave(file.path(pdir, paste0("GO_terms_bins_", binsize, "_", sample1, "_", sample2, "_", status2, ".pdf")), 
	       plot = plotS2, device = cairo_pdf, width = 8, height = 6)
}

