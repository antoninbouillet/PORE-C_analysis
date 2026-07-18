#!/bin/bash

pDir=/home/anton/Bureau/PORE-C_repo/plots

for plot in $pDir/*.pdf; do
	pdfcrop --margins 25 "$plot" "$plot"
done

convert -density 35 "$pDir"/fig1_500kb.pdf "$pDir"/fig1_500kb.tiff
convert -density 68 "$pDir"/fig2_500kb.pdf "$pDir"/fig2_500kb.tiff
convert -density 35 "$pDir"/figS3_500kb.pdf "$pDir"/figS3_500kb.tiff
