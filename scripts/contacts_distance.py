#!/home/anton/venv/cooler/bin/python3
import warnings
import os
from os.path import join
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
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
samples = ["alpha", "beta", "STM"]

binsizes = {"1Mb": 1e6, "500kb": 5e5, "250kb": 2.5e5, "100kb": 1e5}
xlim = (0, 10e7)
ylim = (1e-05, 7e-04)

vmin = 1e2
vmax = 1e6


chromsizes = pd.read_csv(
    join(baseDir, "data/regions/chromsizes.bed"), sep="\t", header=None
)
chr_viewframe = pd.read_csv(join(baseDir, "data/regions/chr_viewframe.bed"), sep="\t")
bioframe.core.checks.is_viewframe(chr_viewframe)


def plotCvd(sample, binsize):

    pDir = join(baseDir, "plots/contacts_distance/%s" % sample)
    minDist = binsizes[binsize]

    if sample == "alpha":
        label = "α"
    elif sample == "beta":
        label = "β"
    else:
        label = sample

    clr = cooler.Cooler(join(cmDir, sample, "%s_%s.cool" % (sample, binsize)))

    # smoothed interactions

    cvd_smooth_agg = cooltools.expected_cis(
        clr=clr,
        view_df=chr_viewframe,
        smooth=True,
        aggregate_smoothed=True,
        smooth_sigma=0.1,
        nproc=num_cpus,
    )

    cvd_smooth_agg["balanced.avg.smoothed"].loc[cvd_smooth_agg["dist"] < 2] = np.nan
    f, ax = plt.subplots(1, 1)
    legend_lines = []
    for region in chr_viewframe["name"]:
        line = ax.loglog(
            cvd_smooth_agg["dist_bp"].loc[cvd_smooth_agg["region1"] == region],
            cvd_smooth_agg["balanced.avg.smoothed"].loc[
                cvd_smooth_agg["region1"] == region
            ],
            label=region,
        )[0]
        legend_lines.append(
            Line2D([0], [0], color=line.get_color(), lw=2, label=region)
        )
        ax.set(
            xlabel="separation, bp", ylabel="IC contact frequency",
            xlim=(2*minDist, 8e7)
            )
        plt.title("%s %s , smoothed" % (label, binsize))
        ax.set_aspect(1.0)
        ax.grid(lw=0.5)
    ax.legend(handles=legend_lines, title="")
    smoothFname = "%s_%s_cvd.pdf" % (sample, binsize)
    plt.savefig(join(pDir, smoothFname), dpi=250, format="pdf")

    cvd_merged = cvd_smooth_agg.drop_duplicates(subset=["dist"])[
        ["dist_bp", "balanced.avg.smoothed.agg"]
    ]

    # Calculate derivative in log-log space

    der = np.gradient(
        np.log(cvd_merged["balanced.avg.smoothed.agg"]), np.log(cvd_merged["dist_bp"])
    )

    f, axs = plt.subplots(
        figsize=(6.5, 13), nrows=2, gridspec_kw={"height_ratios": [4, 2]}, sharex=True
    )
    ax = axs[0]
    ax.loglog(cvd_merged["dist_bp"], cvd_merged["balanced.avg.smoothed.agg"], "-")

    ax.set(ylabel="IC contact frequency", xlim=(3*minDist, 1e8))
    ax.set_aspect(1.0)
    ax.grid(lw=0.5)

    ax = axs[1]
    ax.semilogx(cvd_merged["dist_bp"], der, alpha=0.5)
    ax.set(xlabel="separation, bp", ylabel="slope")
    ax.grid(lw=0.5)
    dFname = "%s_cdv_derivative_%s.pdf" % (sample, binsize)
    plt.savefig(join(pDir, dFname), dpi=250, format="pdf")

    # proportions de contact trans-chromosomiques

    # average contacts, in this case between pairs of chromosomal arms:
    ac = cooltools.expected_trans(
        clr, view_df=None, nproc=num_cpus  # full chromosomes as the view
    )
    # pivot a table to generate a pair-wise average interaction heatmap:
    acp = ac.pivot_table(
        values="balanced.avg", index="region1", columns="region2", observed=True
    ).reindex(index=clr.chromnames, columns=clr.chromnames)
    fs = 14
    f, axs = plt.subplots(
        figsize=(6.0, 5.5),
        ncols=2,
        gridspec_kw={"width_ratios": [20, 1], "wspace": 0.1},
    )
    # assign heatmap and colobar axis:
    ax, cax = axs
    # draw a heatmap, using log-scale for interaction freq-s:
    acpm = ax.imshow(
        acp, cmap="YlOrRd", norm=colors.LogNorm(),  # vmin=vmin, vmax=vmax),
        aspect=1.0
    )
    # assign ticks and labels (ordered names of chromosome arms):
    ax.set(
        xticks=range(len(clr.chromnames)),
        yticks=range(len(clr.chromnames)),
        title=sample,
    )
    ax.set_xticklabels(
        chr_viewframe.name, rotation=30, horizontalalignment="right", fontsize=fs
    )
    ax.set_yticklabels(chr_viewframe.name, fontsize=fs)
    # draw a colorbar of interaction frequencies for the heatmap:
    f.colorbar(acpm, cax=cax, label="IC contact frequency")

    # draw a grid around values:
    ax.set_xticks([x - 0.5 for x in range(1, len(clr.chromnames))], minor=True)
    ax.set_yticks([y - 0.5 for y in range(1, len(clr.chromnames))], minor=True)
    ax.grid(which="minor")

    matrixFname = "%s_trans_matrix_%s.pdf" % (sample, binsize)
    plt.savefig(join(pDir, matrixFname), dpi=250, bbox_inches="tight", format="pdf")


def main():
    for sample in samples:
        for binsize in binsizes:
            plotCvd(sample, binsize)


if __name__ == "__main__":
    main()
