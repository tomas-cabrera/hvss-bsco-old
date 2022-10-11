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
params = ["N", "rv", "Z"]

# Specify kgroup
kgroup = None
kgroup = rustics.INFO_BSE_K_GROUPS.iloc[0]

# Toggle weighting by MW GCs
weight_mw = False

# Initialize hist. generator
hg = rustics.HistogramGenerator("mf", [0, np.inf], nprocs=4, kgroup=kgroup)

# Generate histograms
df_models["hists"] = hg.generate_catalog_histograms(df_models, pm_pbar=True)

# Flatten histogram column, and scale down for legibility
cunits = 4
df_models["hists"] = df_models["hists"].apply(lambda x: x[0] / 10 ** cunits)

# Disable weighting by MW GCs
if not weight_mw:
    df_models["weights_mw"] = [1] * df_models.shape[0]

# Iterate over parameters
for pi, p in enumerate(params):
    print("\tparam=%s" % p)

    # Iterate over unique values of the parameter
    for x in np.sort(df_models[p].unique()):
        # Plot histogram
        df_p = df_models[df_models[p] == x]
        plt.hist(
            df_p["hists"],
            bins=rustics.BINS["Nej"] / 10 ** cunits,
            weights=df_p["weights_mw"],
            label=str(x),
            histtype="step",
            lw=2,
            alpha=0.7,
        )

    # Format plot
    plt.xlabel(r"Number of ejections ($\times 10^%d$)" % cunits)
    plt.ylabel("Number of models")
    plt.ylim(rustics.HIST_LIMS["Nej"])
    plt.legend(
        title=rustics.HEADERS_TO_LABELS[p],
        loc="upper right",
        frameon=False,
    )
    plt.tight_layout()

    # Save and close
    fname = "_".join(("nej_hists", p))
    if type(kgroup) != type(None):
        fname = "_".join((fname, kgroup.initials))
    if weight_mw:
        fname = "_".join((fname, "wmw"))
    fname = ".".join((fname, "pdf"))
    plt.savefig("/".join((rustics.PATH_TO_FIGURES, fname)))
    plt.close()
