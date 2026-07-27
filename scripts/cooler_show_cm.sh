#!/bin/bash

source /home/anton/venv/cooler/bin/activate

baseDir=/home/anton/Bureau/PORE-C_repo_Sanger
cmDir="${baseDir}"/data/contact_maps
pDir="${baseDir}"/plots/contact_maps

binsizes=(1Mb 500kb 250kb)
samples=(alpha beta STM)
chromosomes=(NC_088853.1 NC_088854.1 NC_088855.1 NC_088856.1 NC_088857.1 NC_088858.1 NC_088859.1 NC_088860.1 NC_088861.1 NC_088862.1)
chromsizes=(76070991 61469542 61039741 57946171 57274926 56905015 53672946 51133819 50364239 37310742)

for binsize in "${binsizes[@]}"; do

	for sample in "${samples[@]}"; do

		contactMap="${cmDir}"/"$sample"/"${sample}"_"${binsize}".cool

		for i in "${!chromosomes[@]}"; do

			chromosome=${chromosomes[$i]}
			chromsize=${chromsizes[$i]}

			cooler show "$contactMap" \
			"$chromosome":0-"$chromsize" \
			-s log2 \
			-o "$pDir"/"$sample"_"$binsize"_"$chromosome".pdf

		done

		wait

		pdftk "$pDir"/"$sample"_"$binsize"_*.pdf  output "${pDir}"/"${sample}"_"${binsize}".pdf

		for chromosome in "${chromosomes[@]}"; do
			find "${pDir}" -type f -name "*${chromosome}*.pdf" -exec rm {} ";"
		done

	done
done

for file in "$pDir"; do
       	pdfcrop --margins 25 "$file" "$file"
done
