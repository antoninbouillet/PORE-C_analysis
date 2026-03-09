# Scripts directory

## cooler contact maps

cooler documentation : https://cooler.readthedocs.io/en/latest

- at resolutions of 1Mb, 500kb, 250kb and 100kb
- generates a bed file of genomic bins
- generates contact maps for alpha, beta and STM

## fanc analysis

FAN-C documentation : https://fan-c.readthedocs.io/en/latest

- per chromosome, calculates : 
    - A/B compartments and the corresponding eingenvectors
    - insulation scores

- plots the following contact maps :
   - interaction frequency + insulation scores
   - observed / expected contacts
   - compartments and eingenvectors + gene density 

## fanc compare

- compares the interaction frequency and insulation scores between alpha / beta and alpha / STM
    - the data is log2 transformed after comparison
    - contact maps are normalized to the same number of pairs before comparison

- per chromosome / comparison, plots :
    - a contact map with the log2 fold change of the interaction frequency
    - a heatmap with the log2 fold change of the insulation score
    - the regions with significant interaction differences (detected with HiCcomapre)

## HiCcompare analysis

HiCcompare documentation : https://www.bioconductor.org/packages/release/bioc/html/HiCcompare.html

- between alpha / beta and alpha / STM performs comparisons of interaction frequency for each bin pair
- pairs of bins with a low total interaction frequency (sample1 + sample2 = A value) are filtered
- matrices are normalized to account for the distance effect and differences in coverage between samples
- writes a bed file of regions corresponding to bin pairs with significant interaction differences and a bedpe file of all bins

## HiCcompare intersects

- writes a bed file of all the genes in bin pairs with significant interaction differences (detected with HiCcomapre)

## HiCcomapre GO

- GO enrichement analysis of the genes in regions with significant interaction differences :
    - for all genes
    - for genes in bins that have more contacts in sample1 (alpha) or sample2 (beta / STM)

## HiCcompare plots

- per comparison, draws : 
    - a volcano plot (adjusted M value = log2 fold change vs -log10(P-adj))
    - contact maps of the adjusted M value and -log10(p-adj)
    - adjusted interaction frequency vs distance (log / log), per sample
    - for bin pairs with significant interaction differences, comparison of the average distance between pairs with more contacts in sample1 / sample2
