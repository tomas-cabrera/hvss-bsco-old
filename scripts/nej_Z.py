import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import hvss_utils as hvss
print(hvss.__file__)

###############################################################################

plt.style.use(hvss.MPL_STYLE_FILE_PATH)

# Get hvss_utils info 
Z_list = hvss.CATALOG_MODELS.Z.unique()
bins = hvss.NejBINS

df = pd.read_csv("/".join((hvss.DATA_PATH,"catalog","cmc-catalog.final_dyn.dat")))

for zi, Z in enumerate(Z_list):
    dfz = df[df.Z == Z]
    plt.hist(
        dfz["N_ej"],
        bins=bins,
        label="%g" % (Z/0.02),
        histtype="step",
        lw=2,
        alpha=0.7,
    )
plt.xlabel("Number of ejections/model")
plt.ylabel("Number of models")
plt.xscale("log")
plt.legend(
    title=r"$Z (Z_\odot)$",
)
plt.tight_layout()
plt.savefig("".join((__file__.split(".")[0], ".pdf")))
plt.close()
