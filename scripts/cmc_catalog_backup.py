import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mplticker

import rustics

################################################################################

plt.style.use(rustics.PATH_TO_MPLRC)

# Get model info
df_models = rustics.DF_MODELS

# Iterables
columns = ["vout", "dist", "mf"]
params = ["N", "rv", "Z"]

# Specify kgroup
kgroup = None
kgroup = rustics.INFO_BSE_K_GROUPS.iloc[0]

# Toggle weighting by MW GCs
weight_mw = False

# Make subplot mosiac
plot_size = 2.3
mosaic = [["_".join((p, c)) for c in columns] for p in params]
fig, axd = plt.subplot_mosaic(
    mosaic,
    sharey=True,
    figsize=(plot_size * len(columns), plot_size * len(params)),
    gridspec_kw={
        "wspace": 0,
        "hspace": 0,
    },
)

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

        # Choose axis
        ax = axd["_".join((p, c))]

        # Iterate over unique values of the parameter
        for x in np.sort(df_models[p].unique()):
            # Plot histogram
            mask = df_models[p] == x
            nclusters = mask.sum()
            weights = df_models[mask]["hists"].sum()
            ax.hist(
                bins[:-1],
                bins=bins,
                weights=weights / rustics.N_REALZ / nclusters,
                label=str(x),
                histtype="step",
                lw=2,
                alpha=0.7,
                rasterized=True,
            )

        ## Format plot
        # Set log scale, make sure all ticks are present
        if bins[1] - bins[0] != bins[2] - bins[1]:
            ax.set_xscale("log")
        ax.set_yscale("log")
        xlim = ax.get_xlim()
        numticks = np.ceil(np.log10(xlim[1] / xlim[0])) + 2
        # locmaj = mplticker.LogLocator(base=10, numticks=numticks)
        # ax.xaxis.set_major_locator(locmaj)
        locmin = mplticker.LogLocator(
            base=10, subs=np.arange(0, 1, 0.1), numticks=numticks
        )
        ax.xaxis.set_minor_locator(locmin)
        ax.xaxis.set_minor_formatter(mplticker.NullFormatter())

        # Adjust x-axis labels and ticks; add legend to top row
        if p == params[0]:
            ax.tick_params(labeltop=True)
        elif p == params[-1]:
            ax.set_xlabel(rustics.HEADERS_TO_LABELS[c])
        else:
            ax.set_xticklabels([])

        # Adjust y-axis labels and ticks
        if c == columns[0]:
            ax.set_ylabel(r"Counts/$N_{\rm models}$")
            ax.legend(
                title=rustics.HEADERS_TO_LABELS[p],
                loc=rustics.LOC_LEGEND[c],
                frameon=False,
            )
        elif c == columns[-1]:
            ax.tick_params(labelright=True)
        else:
            ax.set_yticklabels([])

# Make it neat
plt.tight_layout()

# Save and close
fname = "cmc_catalog"
if type(kgroup) != type(None):
    fname = "_".join((fname, kgroup.initials))
if weight_mw:
    fname = "_".join((fname, "wmw"))
fname = ".".join((fname, "pdf"))
plt.savefig("/".join((rustics.PATH_TO_FIGURES, fname)))
plt.close()
