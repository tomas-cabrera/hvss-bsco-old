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
columns = ["vout", "dist", "mf"]
kgroups = [kg for _, kg in rustics.INFO_BSE_K_GROUPS.iterrows()]

# Toggle weighting by MW GCs
weight_mw = True

# Make subplot mosaic
plot_size = 2.3
mosaic = [["_".join((kgroup["label"], c)) for c in columns] for kgroup in kgroups]
fig, axd = plt.subplot_mosaic(
    mosaic,
    sharey=True,
    figsize=(plot_size * len(columns), 0.5 * plot_size * len(kgroups)),
    gridspec_kw={
        "wspace": 0,
        "hspace": 0,
    },
)

# Iterate over all kgroups
for kgroup in kgroups:
    print("kgroup=%s" % kgroup["label"])
    # Iterate over columns (each should have an entry in rustics.BINS_CC)
    for ci, c in enumerate(columns):
        print("column=%s" % c)

        # Initialize hist. generator
        bins = rustics.BINS_CC[c]
        hg = rustics.HistogramGenerator(c, bins, nprocs=4, kgroup=kgroup, cc=True)

        # Generate histograms
        df_models["hists"] = hg.generate_catalog_histograms(df_models, pm_pbar=True)

        # Weight by MW GCs
        if weight_mw:
            df_models["hists"] *= df_models["weights_mw"]

        # Choose axis
        ax = axd["_".join((kgroup["label"], c))]

        # Iterate over unique values of the parameter
        weights = df_models["hists"].sum()
        for cci, cc_row in rustics.INFO_CC_GROUPS.iterrows():
            # Plot histogram
            ax.hist(
                bins[:-1],
                bins=bins,
                weights=weights[cci] / rustics.N_REALZ,
                label=cc_row["plural"],
                histtype="step",
                lw=2,
                alpha=0.7,
            )

        ## Format plot
        # Set log scale, make sure all ticks are present
        if bins[1] - bins[0] != bins[2] - bins[1]:
            ax.set_xscale("log")
        ax.set_yscale("log")

        # Adjust x-axis labels and ticks; add legend to top row
        kgl = kgroup["label"]
        if kgl == kgroups[0]["label"]:
            ax.tick_params(labeltop=True)
        elif kgl == kgroups[-1]["label"]:
            ax.set_xlabel(rustics.HEADERS_TO_LABELS[c])
        else:
            ax.set_xticklabels([])

        # Adjust y-axis labels and ticks
        if c == columns[0]:
            ax.set_ylabel("Counts (%s)" % kgroup.initials)
            # TODO: move legend outside of plot area
            ax.legend(
                loc="lower center",
                frameon=False,
            )
        elif c == columns[-1]:
            ax.tick_params(labelright=True)
        else:
            ax.set_yticklabels([])

# Make it neat
plt.tight_layout()

# Save and close
fname = "cmc_catalog_cc"
if weight_mw:
    fname = "_".join((fname, "wmw"))
fname = ".".join((fname, "pdf"))
plt.savefig("/".join((rustics.PATH_TO_FIGURES, fname)))
plt.close()
