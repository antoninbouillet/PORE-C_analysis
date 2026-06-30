#!/usr/bin/bash

source /home/anton/venv/chess/bin/activate

binsizes=(1Mb 500kb 250kb 100kb)
vminCp=-0.1
vmaxCp=0.1

vminOe=-2
vmaxOe=2

samples=(alpha beta STM)
genome=/home/anton/Bureau/MG_genome/GCF_902806645.1/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fa
baseDir=/home/anton/Bureau/PORE-C_repo
regionsDir=${baseDir}/data/regions

chromosomes=(NC_047559.1 NC_047560.1 NC_047561.1 NC_047562.1 NC_047563.1 NC_047564.1 NC_047565.1 NC_047566.1 NC_047567.1 NC_047568.1)
chromsizes=(55785328 73222313 58319100 53127865 73550375 60151564 62107823 58462999 37089910 57541580)

for binsize in "${binsizes[@]}"; do

	for sample in "${samples[@]}"; do
		
		cmDir=${baseDir}/data/contact_maps/${sample}
		pdir=${baseDir}/plots/fanc_analysis/${sample}
		contact_map=${cmDir}/${sample}_${binsize}.cool
		compartmentFile=${baseDir}/data/matrices/${sample}_${binsize}_compartments.ab

		# calculate compartments
		fanc compartments -g ${genome} -f --recalculate ${contact_map} ${compartmentFile}
		# calculate AB eingenvector
		fanc compartments -f -v ${baseDir}/data/matrices/${sample}_${binsize}_ev.txt ${compartmentFile}
		# calculate insulation scores
		fanc insulation -i ${contact_map} ${baseDir}/data/matrices/${sample}_${binsize}.insulation

		# plots

		for i in "${!chromosomes[@]}"; do

		    chromosome=${chromosomes[$i]}
		    chromsize=${chromsizes[$i]}

		    echo $chromosome

		    fancplot --width 5 -o ${pdir}/${sample}_${binsize}_${chromosome}_oe.pdf \
			    $chromosome:0-${chromsize} -p triangular --title ${chromosome} -e ${contact_map} -vmin -2 -vmax 2

		    fancplot --width 8 -o ${pdir}/${sample}_${binsize}_${chromosome}_compartments.pdf \
			    $chromosome -p square --title ${chromosome} ${compartmentFile} \
			    -vmin ${vminCp} -vmax ${vmaxCp} -c RdBu_r \
			    -p line -f ${baseDir}/data/matrices/${sample}_${binsize}_ev.txt \
			    -p layer ${regionsDir}/genes.bed

		    fancplot --width 6 -o ${pdir}/${sample}_${binsize}_${chromosome}_insulation.pdf \
			    $chromosome -p triangular --title ${chromosome} ${contact_map} -l \
			    -p bar --title "signif (α vs β)" ${regionsDir}/all_signif_regions_alpha_beta_${binsize}.bed \
			    -p bar --title "signif (α vs α+STM)" ${regionsDir}/all_signif_regions_alpha_STM_${binsize}.bed \
			    # -p layer ${regionsDir}/genes.bed
		done

		pdftk ${pdir}/${sample}_${binsize}_*_oe.pdf output ${pdir}/${sample}_${binsize}_oe.pdf
		pdftk ${pdir}/${sample}_${binsize}_*_compartments.pdf output ${pdir}/${sample}_${binsize}_compartments.pdf
		pdftk ${pdir}/${sample}_${binsize}_*_insulation.pdf output ${pdir}/${sample}_${binsize}_insulation.pdf

		for chromosome in "${chromosomes[@]}"; do
			find ${pdir} -type f -name "*${chromosome}*.pdf" -exec rm {} ";"
		done
	done
done


			    #-p scores ${baseDir}/data/matrices/${sample}_${binsize}.insulation \
