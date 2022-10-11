import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import rustics

################################################################################

plt.style.use(rustics.PATH_TO_MPLRC)

# Get model info
df_models = rustics.DF_MODELS
#df_models = df_models.iloc[:10,:]

# Toggle MW GC weighting
weight_mw = True

# Initialize stackedbar generator
sbg = rustics.StackedBarGenerator(chunksize=None, nprocs=4, discrim_ejections=True)

# If histograms have been saved:
fname = "df_models_stackedbar.npy"
if os.path.exists(fname):
    temp = np.load(fname)
    df_models["stackedbar"] = [temp[i,:,:,:] for i in range(temp.shape[0])]
else:
# Generate histograms
    temp = sbg.generate_catalog_stackedbars(df_models, pm_pbar=True)
    np.save(fname, np.stack(temp))
    df_models["stackedbar"] = temp 

# Apply MW GC weights mask
if weight_mw:
    df_models["stackedbar"] *= df_models["weights_mw"]

# Build up stacked bars, one type_f at a time
bottoms = [0] * (rustics.INFO_TYPES_I.shape[0] - 1)
cunits = 7
sbs = df_models["stackedbar"].sum() / 10 ** cunits
if sbg.discrim_ejections:
    ti_totals = sbs.sum(axis=(1, 2))
else:
    ti_totals = sbs.sum(axis=1)
print("ti_totals:", ti_totals)
kwargs_ej = {
    #"edgecolor": "k",
}
kwargs_noej = {
    "hatch": "//",
    "edgecolor": "w",
    "alpha": 0.5,
}
legend_artists = []
legend_labels = []
for tf, tf_row in rustics.INFO_TYPES_F.iterrows():
    if sbg.discrim_ejections:
        for ei in [0, 1]:
            y = sbs[1:, tf, ei]
            print("(tf,ei)=(%d,%d): " % (tf, ei), y)
            if ei:
                bars = plt.bar(
                    rustics.INFO_TYPES_I.loc[1:, "label"],
                    y,
                    bottom=bottoms,
                    label=tf_row.label,
                    color=tf_row.color,
                    **kwargs_ej,
                )
            else:
                bars = plt.bar(
                    rustics.INFO_TYPES_I.loc[1:, "label"],
                    y,
                    bottom=bottoms,
                    label=tf_row.label,
                    color=tf_row.color,
                    **kwargs_noej,
                )
            for bi, b in enumerate(bars):
                # Add 1 to bi, so it can index into TYPES_I correctly
                bi += 1
                if b.get_height() / ti_totals.max() > 0.08:
                    # Make percent label.  text.usetex is enabled in Carl's mplrc file, so the percent sign must be both Python- and LaTeX-escaped
                    pctstr = "%.1f\%%" % (b.get_height() / ti_totals[bi] * 100)
                    plt.text(
                        b.get_x() + b.get_width() / 2.0,
                        b.get_y() + b.get_height() / 2.0,
                        pctstr,
                        ha="center",
                        va="center",
                    )
            bottoms += y
            if ei:
                legend_artists.append(bars)
                legend_labels.append(tf_row.label)
    else:
        y = sbs[1:, tf]
        print("tf=%d: " % tf, y)
        bars = plt.bar(
            rustics.INFO_TYPES_I.loc[1:, "label"],
            y,
            bottom=bottoms,
            label=tf_row.label,
            color=tf_row.color,
        )
        for bi, b in enumerate(bars):
            # Add 1 to bi, so it can index into TYPES_I correctly
            bi += 1
            if b.get_height() / ti_totals.max() > 0.08:
                # Make percent label.  text.usetex is enabled in Carl's mplrc file, so the percent sign must be both Python- and LaTeX-escaped
                pctstr = "%.1f\%%" % (b.get_height() / ti_totals[bi] * 100)
                plt.text(
                    b.get_x() + b.get_width() / 2.0,
                    b.get_y() + b.get_height() / 2.0,
                    pctstr,
                    ha="center",
                    va="center",
                )
        bottoms += y
plt.xlabel("Encounter type")
plt.ylabel(r"Counts ($10^{%d}$)" % cunits)
if sbg.discrim_ejections:
    de_color = "gray"
    legend_artists.append(mpatches.Patch(facecolor=de_color, **kwargs_ej))
    legend_labels.append("Ejection")
    legend_artists.append(mpatches.Patch(facecolor=de_color, **kwargs_noej))
    legend_labels.append("No ejection")
    plt.legend(
        handles=legend_artists,
        labels=legend_labels,
    )
else:
    plt.legend(
        title="Result",
    )
plt.tight_layout()

# Save and close
fname = "encio_stackedbar_catalog"
if weight_mw:
    fname = "_".join((fname, "wmw"))
if sbg.discrim_ejections:
    fname = "_".join((fname, "dej"))
fname = ".".join((fname, "pdf"))
plt.savefig("/".join((rustics.PATH_TO_FIGURES, fname)))
plt.close()
