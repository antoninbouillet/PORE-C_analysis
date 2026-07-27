#!/bin/bash

baseDir=/home/anton/Bureau/PORE-C_repo_Sanger/scripts

Rscript "$baseDir"/HiCcompare_analysis.R
wait
echo "----- running differential analysis ------"
"$baseDir"/HiCcompare_insersect.sh
wait
echo "----- generating plots ------"
Rscript "$baseDir"/HiCcompare_plots.R
wait
echo "----- running GO analysis ------"
Rscript "$baseDir"/HiCcompare_GO.R
wait
echo "----- arranging figures ------"
Rscript "$baseDir"/fig2_arrange.R
wait
echo "----- formatting figures ------"
"$baseDir"/format_fig.sh
