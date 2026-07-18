#!/home/anton/venv/cooler/bin/python3
import warnings
import os
from os.path import join
import subprocess
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import NullLocator
from multiprocessing import Pool
import numpy as np
import pandas as pd
import bioframe
import cooler
import cooltools
from packaging import version

if version.parse(cooltools.__version__) < version.parse("0.5.2"):
    raise AssertionError(
        "tutorial relies on cooltools version 0.5.2 or higher,"
        + "please check your cooltools version and update to the latest"
    )
warnings.filterwarnings("ignore")

# count cpus
num_cpus = os.getenv("SLURM_CPUS_PER_TASK")
if not num_cpus:
    num_cpus = os.cpu_count()
num_cpus = int(num_cpus)

# import des fichiers
baseDir = "/home/anton/Bureau/PORE-C_repo"
cmDir = join(baseDir, "data/contact_maps/")

# old colors
# alpha : #5948d9
# beta : #48D959
# STM : #D95948

# new colors
# alpha : #6BB0ED
# beta : #193D8A
# STM : #EDD472

samplesToColor = {"alpha": "#6BB0ED", "beta": "#193D8A", "STM": "#EDD472"}
binsizes = {"1Mb": 1e6, "500kb": 5e5, "250kb": 2.5e5, "100kb": 1e5}

colorLegend = [
    mpatches.Patch(facecolor=col, alpha=0.7, label=name)
    for name, col in [("α", samplesToColor["alpha"]), ("β", samplesToColor["beta"]), ("α+STM", samplesToColor["STM"])]
]

nameToChr = {"chr1": "NC_047559.1", "chr2": "NC_047560.1",  "chr3": "NC_047561.1", "chr4": "NC_047562.1", "chr5": "NC_047563.1",
             "chr6": "NC_047564.1", "chr7": "NC_047565.1",  "chr8": "NC_047566.1", "chr9": "NC_047567.1", "chr10": "NC_047568.1"}

pDir = join(baseDir, "plots/contacts_distance")

chromsizes = pd.read_csv(
    join(baseDir, "data/regions/chromsizes.bed"), sep="\t", header=None
)
chr_viewframe = pd.read_csv(join(baseDir, "data/regions/chr_viewframe.bed"), sep="\t")
bioframe.core.checks.is_viewframe(chr_viewframe)

# number of bins per chromosome (1 boxplot per bin)
num_bins = 10

# number of contacts to sample from each condition

fixedCount = 1.5e5

# bin distances and compare the interaction frequency between contitions, per chromosome

# run 3 times to get an idea of the variability due to sampling
for i in range(1, 4):
    for binsize in binsizes.keys():

        # first, dowsample to the same number of contacts per sample (cis + trans)
        # note : can also be sampled to the same number of cis contacts (cis_count)

        for sample in samplesToColor.keys():

            clr = cooler.Cooler(join(cmDir, sample, "%s_%s.cool" % (sample, binsize)))

            p = Pool(num_cpus)
            normClrPath = join(cmDir, sample, "%s_%s_fixed.cool" % (sample, binsize))
            cooltools.sample(clr, out_clr_path=normClrPath, cis_count=fixedCount, exact=True, nproc=num_cpus)
            clr_fixed = cooler.Cooler(normClrPath)
            cooler.balance_cooler(clr_fixed, map=p.map, store=True, min_nnz=0)
            p.close()
            p.terminate()

        for region in chr_viewframe["name"]:

            fig, axs = plt.subplots(ncols=num_bins, figsize=(10, 6))

            chrName = nameToChr[region]

            for sampleIdx, sample in enumerate(samplesToColor.keys()):

                clr = cooler.Cooler(join(cmDir, sample, "%s_%s_fixed.cool" % (sample, binsize)))

                cvd_smooth_agg = cooltools.expected_cis(
                    clr=clr,
                    view_df=chr_viewframe,
                    smooth=True,
                    aggregate_smoothed=True,
                    smooth_sigma=0.1,
                    nproc=num_cpus,
                )

                cvd_smooth_agg.loc[cvd_smooth_agg["dist"] < 2, "balanced.avg.smoothed"] = np.nan
                mask = cvd_smooth_agg["region1"] == region
                dist = cvd_smooth_agg.loc[mask, "dist_bp"].values
                val = cvd_smooth_agg.loc[mask, "balanced.avg.smoothed"].values
                count = cvd_smooth_agg['count.sum'].sum().astype(int)

                # Remove NaNs/zero
                valid = ~np.isnan(dist) & ~np.isnan(val) & (val > 0)
                dist = dist[valid]
                val = val[valid]

                min_dist = dist.min()
                max_dist = dist.max()

                if sampleIdx == 0:
                    min_val = val.min()
                    max_val = val.max()
                    max_count = count.max()
                else:
                    min_val = min(min_val, val.min())
                    max_val = max(max_val, val.max())
                    max_count = max(max_count, count.max())

                # Create log-spaced edges
                bin_edges = np.logspace(np.log10(min_dist), np.log10(max_dist), num_bins + 1)

                # Centers (geometric mean of edges)
                bin_centers = np.sqrt(bin_edges[:-1] * bin_edges[1:])
                widths = 0.4 * bin_centers

                # Assign to bins
                bin_idx = np.digitize(dist, bin_edges) - 1
                bin_idx = np.clip(bin_idx, 0, len(bin_centers) - 1)

                # Group values
                groups = [val[bin_idx == i] for i in range(len(bin_centers))]
                groups = [g for g in groups if len(g) > 0]
                positions = bin_centers[: len(groups)]

                # plot all bins whithin the same figure
                """
                bp = axs[b].boxplot(
                    groups,
                    positions=positions,
                    widths=widths,
                    patch_artist=True,
                    showfliers=False,
                    )
                """

                # make subplots for each bin
                for b in range(num_bins):

                    shade = 0.98 - 0.075 * (b % 2)

                    filVals = val[bin_idx == b]
                    bp = axs[b].boxplot(
                        filVals,
                        positions=[sampleIdx],
                        widths=[0.5],
                        patch_artist=True,
                        showfliers=False,
                        )

                    for median in bp['medians']:
                        median.set_color('black')

                    axs[b].set_yscale("log")
                    axs[b].set_facecolor(str(shade))

                    if b > 0:
                        axs[b].yaxis.set_major_locator(NullLocator())
                        axs[b].yaxis.set_minor_locator(NullLocator())
                    else:
                        axs[b].set_ylabel("Contact Frequency")

                    if b == abs(num_bins / 2) - 1:
                        # axs[b].set_title(f"{binsize} - {chrName}, ({num_bins} bins) - {int(fixedCount * 0.001)}k cis contacts per sample")
                        axs[b].set_title(chrName)

                    axs[b].set_xticks([])
                    axs[b].set_xlabel(f" bin {b}")

                    for patch in bp["boxes"]:
                        patch.set_facecolor(samplesToColor[sample])
                        patch.set_alpha(0.7)

                for ax in axs:
                    ax.set_xlim(-0.5, 2.5)
                    ax.set_ylim(min_val * 0.9, max_val * 1.1)

                print(f"Maximum number of reads in {sample} : {max_count}")

            plt.tight_layout()

            fig.legend(
                handles=colorLegend,
                loc="center right",
                bbox_to_anchor=(0.99, 0.85),
                frameon=True,
                fontsize='small'
            )

            pFname = "cvd_boxplot_%s_%s_run%s.pdf" % (binsize, chrName, i)
            plt.savefig(join(pDir, pFname), dpi=250, bbox_inches="tight", format="pdf")

        bash = "/bin/bash"
        mergePdf = "pdftk %s/cvd_boxplot_%s_NC_*_run%s.pdf output %s/cvd_boxplot_%s_run%s.pdf" % (pDir, binsize, i, pDir, binsize, i)
        subprocess.call(mergePdf, shell=True, executable=bash)

        rm = """find %s -maxdepth 1 -type f -name "*NC_*" -exec rm {} ";" """ % (pDir)
        subprocess.call(rm, shell=True, executable=bash)
