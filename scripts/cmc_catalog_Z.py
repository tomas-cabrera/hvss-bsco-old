import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib import cm
import seaborn as sns
from scipy import stats
from matplotlib.ticker import LogLocator 

import hvss_utils

################################################################################

plt.style.use(hvss_utils.MPL_STYLE_FILE_PATH)

S_SCAT = 0.5

# Get CATALOG_MODELS (note: different from other files!)
catalog_models = hvss_utils.CATALOG_MODELS
N_list = catalog_models.N.unique()
rv_list = catalog_models.rv.unique()
Z_list = catalog_models.Z.unique()

N = 10
vbins = hvss_utils.VBINS
binsingle_types = hvss_utils.BINSINGLE_TYPES_CMC_I
bse_k_groups = hvss_utils.BSE_K_GROUPS
age_groups = hvss_utils.AGE_GROUPS

# File name
data_file_path = hvss_utils.DATA_PATH + "catalog/{}_N-{}".format(
    hvss_utils.FB_OUTPUT_FILENAME, N
)

# read saved data file
try:
    catalog = hvss_utils.load_flattened(
        data_file_path + "_wbinned.txt",
        hvss_utils.CATALOG_SHAPE_Z,
    )
except Exception as e:
    print("\tError in loading full catalog data:")
    print("\t\t", e)


# N histograms
print("\tGenerating N histograms...")
for gi, bse_k_group in enumerate(bse_k_groups):
    # Counter to skip plot formatting if no data are plotted
    plot_exists = 0
    for ni, N_cluster in enumerate(N_list):
        filtered = catalog[ni, :, :, :, gi, :, :].sum(
            axis=(
                0,
                1,
                2,
                3,
            )
        )
        if filtered.sum() > 0:
            plot_exists = 1
            plt.hist(
                vbins[1:],
                bins=vbins,
                weights=filtered / N,
                label="%d" % N_cluster,
                histtype="step",
                lw=2,
                alpha=0.7,
            )
    if plot_exists:
        plt.xlabel(hvss_utils.HEADERS_TO_LABELS["vout"])
        plt.ylabel("Counts")
        plt.ylim(hvss_utils.CATALOG_HIST_VLIM)
        plt.xscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])
        plt.yscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])
        plt.xticks(ticks=[10. ** x for x in np.arange(-2, 5)])
        plt.gca().xaxis.set_minor_locator(LogLocator(subs="all", numticks=100))
        legend = plt.legend(
            title=r"$N [\times 10^5]$",
            loc="upper left",
            frameon=False,
        )
        # for handle in legend.legendHandles:
        #    handle.set_alpha(1)
        plt.savefig(hvss_utils.FIGS_PATH + "cmc_catalog_Z_vouthist_N_{}.pdf".format(bse_k_group.short_name))
        plt.close()
    else:
        print("\t\tNo data found for {}; continuing...".format(bse_k_group.name))
        continue

# rv histograms
print("\tGenerating rv histograms...")
for gi, bse_k_group in enumerate(bse_k_groups):
    # Counter to skip plot formatting if no data are plotted
    plot_exists = 0
    for ri, rv in enumerate(rv_list):
        filtered = catalog[:, ri, :, :, gi, :, :].sum(
            axis=(
                0,
                1,
                2,
                3,
            )
        )
        if filtered.sum() > 0:
            plot_exists = 1
            plt.hist(
                vbins[1:],
                bins=vbins,
                weights=filtered / N,
                label=rv,
                histtype="step",
                lw=2,
                alpha=0.7,
            )
    if plot_exists:
        plt.xlabel(hvss_utils.HEADERS_TO_LABELS["vout"])
        plt.ylabel("Counts")
        plt.ylim(hvss_utils.CATALOG_HIST_VLIM)
        plt.xscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])
        plt.yscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])
        plt.xticks(ticks=[10. ** x for x in np.arange(-2, 5)])
        plt.gca().xaxis.set_minor_locator(LogLocator(subs="all", numticks=100))
        legend = plt.legend(
            title=r"$r_{\rm vir} [pc]$",
            loc="upper left",
            frameon=False,
        )
        # for handle in legend.legendHandles:
        #    handle.set_alpha(1)
        plt.savefig(hvss_utils.FIGS_PATH + "cmc_catalog_Z_vouthist_rv_{}.pdf".format(bse_k_group.short_name))
        plt.close()
    else:
        print("\t\tNo data found for {}; continuing...".format(bse_k_group.name))
        continue

# Z histograms
print("\tGenerating Z histograms...")
for gi, bse_k_group in enumerate(bse_k_groups):
    # Counter to skip plot formatting if no data are plotted
    plot_exists = 0
    for zi, Z in enumerate(Z_list):
        filtered = catalog[:, :, zi, :, gi, :, :].sum(
            axis=(
                0,
                1,
                2,
                3,
            )
        )
        if filtered.sum() > 0:
            plot_exists = 1
            plt.hist(
                vbins[1:],
                bins=vbins,
                weights=filtered / N,
                label=Z,
                histtype="step",
                lw=2,
                alpha=0.7,
            )
    if plot_exists:
        plt.xlabel(hvss_utils.HEADERS_TO_LABELS["vout"])
        plt.ylabel("Counts")
        plt.ylim(hvss_utils.CATALOG_HIST_VLIM)
        plt.xscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])
        plt.yscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])
        plt.xticks(ticks=[10. ** x for x in np.arange(-2, 5)])
        plt.gca().xaxis.set_minor_locator(LogLocator(subs="all", numticks=100))
        legend = plt.legend(
            title=r"$Z$",
            loc="upper left",
            frameon=False,
        )
        # for handle in legend.legendHandles:
        #    handle.set_alpha(1)
        plt.savefig(hvss_utils.FIGS_PATH + "cmc_catalog_Z_vouthist_Z_{}.pdf".format(bse_k_group.short_name))
        plt.close()
    else:
        print("\t\tNo data found for {}; continuing...".format(bse_k_group.name))
        continue

# age histograms
print("\tGenerating age histograms...")
for gi, bse_k_group in enumerate(bse_k_groups):
    # Counter to skip plot formatting if no data are plotted
    plot_exists = 0
    for agi, age in enumerate(age_groups):
        filtered = catalog[:, :, :, :, gi, agi, :].sum(
            axis=(
                0,
                1,
                2,
                3, 
            )
        )
        if filtered.sum() > 0:
            plot_exists = 1
            plt.hist(
                vbins[1:],
                bins=vbins,
                weights=filtered / N,
                label=age.label,
                color=age.color,
                histtype="step",
                lw=2,
                alpha=0.7,
            )
    if plot_exists:
        plt.xlabel(hvss_utils.HEADERS_TO_LABELS["vout"])
        plt.ylabel("Counts")
        plt.ylim(hvss_utils.CATALOG_HIST_VLIM)
        plt.xscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])
        plt.yscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])
        plt.xticks(ticks=[10. ** x for x in np.arange(-2, 5)])
        plt.gca().xaxis.set_minor_locator(LogLocator(subs="all", numticks=100))
        legend = plt.legend(
            title="Core collapse status",
            loc="upper left",
            frameon=False,
        )
        # for handle in legend.legendHandles:
        #    handle.set_alpha(1)
        plt.savefig(hvss_utils.FIGS_PATH + "cmc_catalog_Z_vouthist_age_{}.pdf".format(bse_k_group.short_name))
        plt.close()
    else:
        print("\t\tNo data found for {}; continuing...".format(bse_k_group.name))
        continue
