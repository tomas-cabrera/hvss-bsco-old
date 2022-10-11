import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator

import rustics

###############################################################################

GRONDIN_RVS = np.array(
    [
        -116.314,
        -189.679,
        -160.858,
        -110.683,
        -73.028,
        -137.845,
        -205.337,
        -191.682,
        20.833,
        -199.733,
        -132.903,
        -138.898,
        -67.637,
    ]
)

GRONDIN_M3RV = -147.20

###############################################################################

plt.style.use(rustics.PATH_TO_MPLRC)

# Get model info
df_models = rustics.DF_MODELS
df_models = df_models[
    (df_models.N == 16)
    & (df_models.rv == 1)
    & (df_models.rg == 8)
    & (df_models.Z == 0.0002)
]

# Iterables
columns = [
    "v_radial",
]

# Specify kgroup
kgroup = None
kgroup = rustics.INFO_BSE_K_GROUPS.iloc[0]

# Iterate over columns (each should have an entry in rustics.BINS)
for ci, c in enumerate(columns):
    print("column=%s" % c)

    bins = rustics.BINS[c]
    df = pd.read_csv(
        rustics.PATH_TO_DATA + "/N1.6e6_rv1_rg8_Z0.0002/output_N-10_ejections.txt"
    )
    df["vout"] = ((df["vfin"] * df["v_crit"]) ** 2 - df["vesc"] ** 2) ** 0.5
    df["v_radial"] = (
        df["vout"] * np.cos(np.pi * np.random.rand(df.shape[0])) + GRONDIN_M3RV
    )
    weights, bins = np.histogram(
        df[(df.kf >= kgroup.lo) & (df.kf < kgroup.hi)]["v_radial"], bins
    )

    # Plot histogram
    plt.hist(
        bins[:-1],
        bins=bins,
        weights=weights / rustics.N_REALZ,
        label="N1.6e6_rv1_rg8_Z0.0002",
        histtype="step",
        lw=2,
        alpha=0.7,
    )
    gweights = [
        8,
        12,
        18,
        40,
        55,
        90,
        133,
        210,
        205,
        210,
        215,
        225,
        270,
        297,
        270,
        195,
        123,
        85,
        60,
        53,
        35,
        27,
        20,
        8,
        9,
        5,
    ]
    gbins = np.linspace(-240, 0, len(gweights) + 1)
    print(len(gweights), len(gbins))
    plt.hist(
        gbins[:-1],
        bins=gbins,
        weights=gweights,
        label=r"Grondin+22 $\texttt{corespray}$",
        histtype="step",
        lw=2,
        alpha=0.7,
    )
    ylim = plt.ylim()
    plt.vlines(
        GRONDIN_RVS,
        *ylim,
        colors="k",
        label="Grondin+22 ETS",
        lw=0.5,
    )
    plt.vlines(
        GRONDIN_M3RV,
        *ylim,
        colors="k",
        label="M3",
        lw=0.5,
        ls="--",
    )

    # Format plot
    plt.xlabel(rustics.HEADERS_TO_LABELS[c])
    plt.ylabel("Counts")
    plt.xlim((-300, 100))
    plt.ylim(ylim[0], ylim[1] + 150)
    if bins[1] - bins[0] != bins[2] - bins[1]:
        plt.xscale("symlog", subs=[2, 3, 4, 5, 6, 7, 8, 9])
        # plt.xticks(ticks=[10. ** int(x) for x in np.arange(*np.log10(rustics.AXES_LIMS["vout"]))])
        plt.gca().xaxis.set_minor_locator(LogLocator(subs="all", numticks=100))
    # plt.yscale("log", subs=[2,3,4,5,6,7,8,9])
    plt.legend(
        frameon=False,
    )
    plt.tight_layout()

    # Save and close
    fname = "_".join(("grondin_ets", c))
    if type(kgroup) != type(None):
        fname = "_".join((fname, kgroup.initials))
    fname = ".".join((fname, "pdf"))
    plt.savefig("/".join((rustics.PATH_TO_FIGURES, fname)))
    plt.close()
