#!/bin/bash

# intersect regions with significant interaction differences with genes

baseDir=/home/anton/Bureau/PORE-C_repo
binsize=500kb

comparisons=(alpha_beta alpha_STM)

for comp in "${comparisons[@]}"; do

	bedtools intersect -wa -wb \
	       	-a ${baseDir}/data/regions/all_signif_regions_${comp}_${binsize}.bed \
		-b ${baseDir}/data/regions/genes.bed \
		> ${baseDir}/data/regions/all_signif_genes_${comp}_${binsize}.bed
done

