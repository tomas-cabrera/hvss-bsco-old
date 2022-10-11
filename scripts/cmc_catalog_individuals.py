import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator

import rustics

################################################################################

plt.style.use(rustics.PATH_TO_MPLRC)

# Get model info
df_models = rustics.DF_MODELS

# Iterables
columns = ["vout", "mf", "dist"]
params = ["N", "rv", "Z"]

# Specify kgroup
kgroup = None
kgroup = rustics.INFO_BSE_K_GROUPS.iloc[0]

# Toggle weighting by MW GCs
weight_mw = False

# Iterate over columns (each should have an entry in rustics.BINS)
for ci, c in enumerate(columns):
    print("column=%s" % c)

    # Initialize hist. generator
    bins = rustics.BINS[c]
    hg = rustics.HistogramGenerator(c, bins, nprocs=4, kgroup=kgroup)

    # Generate histograms
    df_models["hists"] = hg.generate_catalog_histograms(df_models, pm_pbar=True)

    # Weight by MW GCs
    if weight_mw:
        df_models["hists"] *= df_models["weights_mw"]

    # Iterate over parameters
    for pi, p in enumerate(params):
        print("\tparam=%s" % p)

        # Iterate over unique values of the parameter
        for x in np.sort(df_models[p].unique()):
            # Plot histogram
            mask = df_models[p] == x
            nclusters = mask.sum()
            weights = df_models[mask]["hists"].sum()
            plt.hist(
                bins[:-1],
                bins=bins,
                weights=weights / rustics.N_REALZ / nclusters,
                label=str(x),
                histtype="step",
                lw=2,
                alpha=0.7,
            )

        # Format plot
        plt.xlabel(rustics.HEADERS_TO_LABELS[c])
        plt.ylabel(r"Counts/$N_{\rm models}$")
        plt.ylim(rustics.HIST_LIMS[c])
        if bins[1] - bins[0] != bins[2] - bins[1]:
            plt.xscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])
            # plt.xticks(ticks=[10. ** int(x) for x in np.arange(*np.log10(rustics.AXES_LIMS["vout"]))])
            plt.gca().xaxis.set_minor_locator(LogLocator(subs="all", numticks=100))
        plt.yscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])
        plt.legend(
            title=rustics.HEADERS_TO_LABELS[p],
            loc=rustics.LOC_LEGEND[c],
            frameon=False,
        )
        plt.tight_layout()

        # Save and close
        fname = "_".join(("cmc_catalog", c, p))
        if type(kgroup) != type(None):
            fname = "_".join((fname, kgroup.initials))
        if weight_mw:
            fname = "_".join((fname, "wmw"))
        fname = ".".join((fname, "pdf"))
        plt.savefig("/".join((rustics.PATH_TO_FIGURES, fname)))
        plt.close()
