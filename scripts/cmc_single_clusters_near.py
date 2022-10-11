import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib import cm
from matplotlib.legend_handler import HandlerTuple
import seaborn as sns
from scipy import stats
from inspect import stack

import rustics

################################################################################

plt.style.use(rustics.PATH_TO_MPLRC)

# Get cluster names
cmc_cluster_list = [
    "N8e5_rv0.5_rg8_Z0.0002",
    "N8e5_rv2_rg8_Z0.0002",
    "N8e5_rv0.5_rg8_Z0.02",
    "N4e5_rv0.5_rg8_Z0.0002",
]
#cmc_cluster_list = rustics.DF_MODELS[(rustics.DF_MODELS.N == 4) & (rustics.DF_MODELS.rv == 0.5)].fname

for ci, cmc_cluster in enumerate(cmc_cluster_list):
    print("{}, N_REALZ={}".format(cmc_cluster, rustics.N_REALZ))

    # File names and directory creation
    path_to_ejections = "/".join(
        (rustics.PATH_TO_DATA, cmc_cluster, rustics.FILE_EJECTIONS)
    )

    # Read saved data file
    ejdf = rustics.EjectionDf(path_to_ejections)

    # Convert to physical units
    ejdf.convert_from_fewbody()

    # Filter out mergers
    ejdf.df = ejdf.df[ejdf.df["type_f"] > 0]

    # Calculate vout (velocity from cluster)
    ejdf.df["vout"] = ejdf.calc_vout()
    ejdf.df["mx"] = ejdf.df[["m0", "m10", "m11"]].max(axis=1)
    ejdf.df["dist"] = (
        ejdf.df.vout
        * (14000. - ejdf.df.time)
        * (1e6 * rustics.YR_TO_S)
        / (rustics.PC_TO_M)
    )

    # TEMP: Trim to objects still w/i 100 pc of their host clusters
    ejdf.df = ejdf.df[ejdf.df.dist < 0.1]

    # Load times and core radii from initial.dyn.dat
    df_dyn = pd.read_csv(
        "/".join((rustics.PATH_TO_CMC_DATA, cmc_cluster, "initial.dyn.dat")),
        names=["t", "N", "rc_nb"],
        usecols=[0, 3, 24],  # 0 for t, 3 for N, 24 for rc_nb
        skiprows=2,
        delim_whitespace=True,
    )

    # Load intial.conv.sh (for converting TotalTime to Myr and core radius to pcs)
    # Copied from cmctools/cmctoolkit.py
    conv_path = "/".join((rustics.PATH_TO_CMC_DATA, cmc_cluster, "initial.conv.sh"))
    f = open(conv_path, "r")
    convfile = f.read().split("\n")
    f.close()
    timeunitsmyr = float(convfile[19][13:])
    lengthunitparsec = float(convfile[15][17:])
    df_dyn["t"] *= timeunitsmyr
    df_dyn["rc_nb"] *= lengthunitparsec

    # Generate ordered list of type_i's, descending by number of occurences
    # type_is = [x[0] for x in ejdf.df.value_counts("type_i").index]
    if ci == 0:
        type_is = ejdf.df.value_counts("type_i").index
        type_is = rustics.INFO_TYPES_I.iloc[type_is]

    # Some parameters
    ms = 0.25
    minimum = 50

    # Select x-/y-axis labels; stars only (k<10)
    x = "time"
    y = "mf"
    kgroup = rustics.INFO_BSE_K_GROUPS.iloc[0]

    # Set up axes/subplots
    axis_ratio = 3
    fig, axd = plt.subplot_mosaic(
        [["thist", "."], ["scatter", "vouthist"]],
        gridspec_kw={
            "hspace": 0.0,
            "wspace": 0.0,
            "height_ratios": [1, axis_ratio],
            "width_ratios": [axis_ratio, 1],
        },
        figsize=(4, 4),
    )
    axd["scatter"].set_xscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])
    axd["scatter"].set_yscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])
    axd["thist"].set_xscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])
    axd["thist"].set_yscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])
    axd["vouthist"].set_xscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])
    axd["vouthist"].set_yscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])

    # Add r_core axis
    colorr = "xkcd:crimson"
    axr = axd["thist"].twinx()
    axr.set_yscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])
    tsample = [
        int(i)
        for i in np.logspace(
            0,
            np.log10(df_dyn.shape[0] - 1),
            10000,
        )
    ]
    tsample = list(set(tsample))
    tsample.sort()
    axr.plot(
        df_dyn.loc[tsample, "t"] / 1000.0,
        df_dyn.loc[tsample, "rc_nb"],
        color=colorr,
        alpha=0.7,
        lw=0.5,
        rasterized=True,
    )
    axr.set_ylabel(r"$r_{c, NB} [{\rm pc}]$", color=colorr)
    axr.tick_params(
        axis="y",
        colors=colorr,
    )
    axr.set_ylim(
        (0.0011, axr.get_ylim()[1])
    )  # hacky way of getting 0.0 label off the adjacent axes

    legend_artists = []
    legend_labels = []
    nbins = 30
    for ti, ti_row in type_is.iterrows():
        # for ti, ti_row in rustics.INFO_TYPES_I.iterrows():
        print(ti, ti_row)
        filtered = ejdf.df[
            (ejdf.df.kf >= kgroup.lo)
            & (ejdf.df.kf < kgroup.hi)
            & (ejdf.df.type_i == ti)
        ]
        print("\tfiltered.shape: ", filtered.shape)
        (s,) = axd["scatter"].plot(
            filtered[x] / 1000.0,
            filtered[y],
            color=ti_row.color,
            linestyle="None",
            marker="+",
            markersize=ms,  # irrelevant for marker="," (pixels)
            alpha=0.25,
            rasterized=True,
        )
        if ci == 0:
            xlims = axd["scatter"].get_xlim()
            ylims = axd["scatter"].get_ylim()
            ylims = (0.05, ylims[1])
        _, _, h, = axd["thist"].hist(
            filtered[x] / 1000.0,
            # density=True,
            bins=np.logspace(*np.log10(xlims), nbins),
            color=ti_row.color,
            histtype="step",
            lw=2,
            alpha=0.7,
            rasterized=True,
        )
        axd["vouthist"].hist(
            filtered[y],
            # density=True,
            bins=np.logspace(*np.log10(ylims), nbins),
            color=ti_row.color,
            orientation="horizontal",
            histtype="step",
            lw=2,
            alpha=0.7,
            rasterized=True,
        )
        legend_artists.append((s, h[0]))
        legend_labels.append(ti_row.label)
    axd["scatter"].set_xlabel(rustics.HEADERS_TO_LABELS[x])
    axd["scatter"].set_ylabel(rustics.HEADERS_TO_LABELS[y])
    axd["thist"].set_xlabel("")
    axd["thist"].set_ylabel(r"$N$")
    axd["vouthist"].set_xlabel(r"$N$")
    axd["vouthist"].set_ylabel("")
    print(axd["vouthist"].get_xticks())
    axd["vouthist"].set_xticklabels(["0", "0.01"])
    # axd["vouthist"].set_xticklabels([str(l) for l in axd["vouthist"].get_xticks()])
    axd["vouthist"].set_yticklabels([])
    if ci == 0:
        ylimh = axd["thist"].get_ylim()
        xlimh = axd["vouthist"].get_xlim()
        xlimh = (0.0, 0.04)
        ylimh = (1, 1e5)
    axd["scatter"].set_xlim(xlims)
    #axd["scatter"].set_ylim(ylims)
    axd["thist"].set_xlim(xlims)
    axd["thist"].set_ylim(ylimh)
    # axd["vouthist"].set_xlim((0.0002, axd["vouthist"].get_xlim()[1]))
    #axd["vouthist"].set_ylim(ylims)

    # Add legend
    legend = axd["scatter"].legend(
        title=cmc_cluster,
        handles=legend_artists,
        labels=legend_labels,
        handler_map={tuple: HandlerTuple(ndivide=None)},
        markerscale=8.0 / ms,
        frameon=False,
        loc="lower left",
        prop={"size": 6.5 * fig.get_figheight() / axis_ratio},
        ncol=2,
    )
    for handle in legend.legendHandles:
        handle.set_alpha(1)

    # Adjust spacing
    plt.tight_layout()

    # hvss_utils.grid_selfmade(axd["scatter"])
    Nstr, rvstr, rgstr, Zstr = cmc_cluster.split("_")
    lstr = "-".join((y, Nstr, rvstr, rgstr, Zstr))
    plt.savefig(
        "/".join(
            (
                rustics.PATH_TO_FIGURES,
                lstr.join(("cmc_single_clusters_temp_", ".pdf")),
            )
        )
    )
    plt.close()
