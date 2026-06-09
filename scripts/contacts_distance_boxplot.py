#!/home/anton/venv/cooler/bin/python3
import warnings
import os
from os.path import join
import subprocess
import matplotlib.pyplot as plt
from matplotlib import colors
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

# Analyse des cartes de contact : fréquence de contact en fonction de la distance
# Source https://cooltools.readthedocs.io/en/latest/notebooks/contacts_vs_distance.html

# count cpus
num_cpus = os.getenv("SLURM_CPUS_PER_TASK")
if not num_cpus:
    num_cpus = os.cpu_count()
num_cpus = int(num_cpus)

# import des fichiers
baseDir = "/home/anton/Bureau/PORE-C_repo"
cmDir = join(baseDir, "data/contact_maps/")
samplesToColor = {"alpha": "#5948d9", "beta": "#48D959", "STM": "#D95948"}
binsizes = {"1Mb": 1e6, "500kb": 5e5, "250kb": 2.5e5, "100kb": 1e5}

pDir = join(baseDir, "plots/contacts_distance")

chromsizes = pd.read_csv(
    join(baseDir, "data/regions/chromsizes.bed"), sep="\t", header=None
)
chr_viewframe = pd.read_csv(join(baseDir, "data/regions/chr_viewframe.bed"), sep="\t")
bioframe.core.checks.is_viewframe(chr_viewframe)

num_bins = 10
binsize = "500kb"

# bin distances and compare the interaction frequency between contitions, per chromosome

for binsize in binsizes.keys():

    for region in chr_viewframe["name"]:

        fig, ax = plt.subplots(figsize=(10, 6))

        for sample in samplesToColor.keys():

            clr = cooler.Cooler(join(cmDir, sample, "%s_%s.cool" % (sample, binsize)))

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

            # Remove NaNs/zero
            valid = ~np.isnan(dist) & ~np.isnan(val) & (val > 0)
            dist = dist[valid]
            val = val[valid]

            min_dist = dist.min()
            max_dist = dist.max()

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

            bp = ax.boxplot(
                groups,
                positions=positions,
                widths=widths,
                patch_artist=True,
                showfliers=False,
            )

            for patch in bp["boxes"]:
                patch.set_facecolor(samplesToColor[sample])
                patch.set_alpha(0.7)

        ax.set_xscale("log")
        ax.set_yscale("log")

        ax.set_xlabel("Separation, bp")
        ax.set_ylabel("IC contact frequency")
        ax.set_title(f"{binsize} - {region}, ({num_bins} bins)")
        ax.grid(lw=0.5, which="both")

        plt.tight_layout()
        pFname = "cvd_boxplot_%s_%s.pdf" % (binsize, region)
        plt.savefig(join(pDir, pFname), dpi=250, bbox_inches="tight", format="pdf")

    bash = "/bin/bash"
    mergePdf = "pdftk %s/cvd_boxplot_%s_chr*.pdf output %s/cvd_boxplot_%s.pdf" % (pDir, binsize, pDir, binsize)
    subprocess.call(mergePdf, shell=True, executable=bash)

    rm = """find %s -maxdepth 1 -type f -name "*chr*" -exec rm {} ";" """ % (pDir)
    subprocess.call(rm, shell=True, executable=bash)
