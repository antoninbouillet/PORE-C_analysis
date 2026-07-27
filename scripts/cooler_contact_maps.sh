#!/bin/bash

samples=(alpha beta STM)
binNames=(1Mb 500kb 250kb 100kb)
binSizes=(1000000 500000 250000 100000)
baseDir=/home/anton/Bureau/PORE-C_repo_Sanger
rDir=${baseDir}/data/regions
cmDir=${baseDir}/data/contact_maps
pairDir=${baseDir}/data/pairs

source /home/anton/venv/cooler/bin/activate

# doc : https://cooler.readthedocs.io/en/latest/
#
for i in "${!binNames[@]}"; do

	binName=${binNames[$i]}
	binSize=${binSizes[$i]}

	# generate bins
	
	cooler makebins ${rDir}/chromsizes.bed ${binSize} -o ${rDir}/bins_${binName}.bed

	# generate contact maps
	
	for sample in "${samples[@]}"; do

		cooler cload pairs -c1 2 -p1 3 \
		-c2 4 -p2 5 \
		${rDir}/bins_${binName}.bed \
		${pairDir}/${sample}/${sample}.pairs.gz \
		${cmDir}/${sample}/${sample}_${binName}.cool

		wait

		cooler balance ${cmDir}/${sample}/${sample}_${binName}.cool
				
	done
	
done
