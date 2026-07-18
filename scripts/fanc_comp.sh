#!/usr/bin/bash

source /home/anton/venv/chess/bin/activate

binsizes=(1Mb 500kb 250kb)
comp=(beta STM)
vmin=-3
vmax=3

genome=/home/anton/Bureau/MG_genome/GCF_902806645.1/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fa
baseDir=/home/anton/Bureau/PORE-C_repo
cmDir=${baseDir}/data/contact_maps
mDir=${baseDir}/data/matrices
regionsDir=${baseDir}/data/regions

# doc : https://fan-c.readthedocs.io/en/latest/fanc-executable/fanc-analyse-hic/comparisons.html 

sample1=alpha
chromosomes=(NC_047559.1 NC_047560.1 NC_047561.1 NC_047562.1 NC_047563.1 NC_047564.1 NC_047565.1 NC_047566.1 NC_047567.1 NC_047568.1)
chromsizes=(55785328 73222313 58319100 53127865 73550375 60151564 62107823 58462999 37089910 57541580)

for binsize in "${binsizes[@]}"; do
	: '
	for sample2 in ${comp[@]}; do

		pdir=${baseDir}/plots/fanc_comparison/${sample1}_${sample2}

		contact_map1=${cmDir}/${sample1}/${sample1}_${binsize}.cool
		contact_map2=${cmDir}/${sample2}/${sample2}_${binsize}.cool

		compartmentFile1=${mDir}/${sample1}_${binsize}_compartments.ab
		compartmentFile2=${mDir}/${sample2}_${binsize}_compartments.ab

		# comparison of counts matrices are normalized to the same number of pairs before comparison
		fanc compare -Z -l -I ${contact_map1} ${contact_map2} ${mDir}/comp_${binsize}_${sample1}_vs_${sample2}.fanc

		# comparison of observed / expected counts : no normalization
		fanc compare -Z -l -I -e ${contact_map1} ${contact_map2} ${mDir}/comp_oe_${binsize}_${sample1}_vs_${sample2}.fanc

		# compare insulation scores
		fanc compare -Z -l -I \
			${mDir}/${sample1}_${binsize}.insulation \
			${mDir}/${sample2}_${binsize}.insulation \
			${mDir}/${sample1}_vs_${sample2}_${binsize}.insulation


		compFile=${mDir}/comp_${binsize}_${sample1}_vs_${sample2}.fanc
		# compFileOe=${mDir}/comp_oe_${binsize}_${sample1}_vs_${sample2}.fanc

		for i in "${!chromosomes[@]}"; do

			chromosome=${chromosomes[$i]}
			chromsize=${chromsizes[$i]}

			echo $chromosome

			fancplot --width 6 -o ${pdir}/${chromosome}_${binsize}_${sample1}_vs_${sample2}_comp.pdf ${chromosome}:0-${chromsize} \
				-p triangular --title "${chromosome} : log2 fold change (contacts)" ${compFile} -vmin ${vmin} -vmax ${vmax} -c RdBu_r \
				-p scores --title "log2 fold change (insulation score)" ${mDir}/${sample1}_vs_${sample2}_${binsize}.insulation \
				-p bar --title "signif (HiCcompare)" ${regionsDir}/all_signif_regions_${sample1}_${sample2}_${binsize}.bed 

			done

			pdftk ${pdir}/*_${binsize}_${sample1}_vs_${sample2}_comp.pdf output ${pdir}/${binsize}_${sample1}_vs_${sample2}_comp.pdf
			# pdftk ${pdir}/*_${binsize}_${sample1}_vs_${sample2}_comp_oe.pdf output ${pdir}/${binsize}_${sample1}_vs_${sample2}_comp_oe.pdf

		for chromosome in "${chromosomes[@]}"; do
			find ${pdir} -type f -name "*${chromosome}*.pdf" -exec rm {} ";"
		done
	done
	'
	# mirror plots : log2FC alpha/beta vs alpha/STM

	pdir=${baseDir}/plots/fanc_comparison

	comp1=${mDir}/comp_${binsize}_alpha_vs_beta.fanc
	comp2=${mDir}/comp_${binsize}_alpha_vs_STM.fanc

	comp_oe1=${mDir}/comp_oe_${binsize}_alpha_vs_beta.fanc
	comp_oe2=${mDir}/comp_oe_${binsize}_alpha_vs_STM.fanc


	for i in "${!chromosomes[@]}"; do

		chromosome="${chromosomes[$i]}"
		chromsize="${chromsizes[$i]}"

		fancplot -o "${pdir}"/mirror_ct_"${chromosome}"_"${binsize}".pdf \
			"${chromosome}":0-"${chromsize}" \
			-p mirror  -lvmin -4.5 -lvmax 4.5 -uvmin -4.5 -uvmax 4.5 \
			-uc RdBu_r -lc RdBu_r "${comp1}" "${comp2}" \
			--title "log2 fold change (α vs β / α vs α + STM)"

		fancplot -o "${pdir}"/mirror_oe_"${chromosome}"_${binsize}.pdf \
			"${chromosome}":0-"${chromsize}" \
			-p mirror -lvmin -4 -lvmax 4 -uvmin -4 -uvmax 4 \
			-uc RdBu_r -lc RdBu_r "${comp_oe1}" "${comp_oe2}" \
			--title "log2 fold change observed / expected (α vs β / α vs α + STM)"

	done

	pdftk "${pdir}"/mirror_ct_*_"${binsize}".pdf output "${pdir}"/mirror_ct_"${binsize}".pdf
	pdftk "${pdir}"/mirror_oe_*_"${binsize}".pdf output "${pdir}"/mirror_oe_"${binsize}".pdf

	for chromosome in "${chromosomes[@]}"; do
		find ${pdir} -type f -name "*${chromosome}_${binsize}.pdf" -exec rm {} ";"
	done

done
