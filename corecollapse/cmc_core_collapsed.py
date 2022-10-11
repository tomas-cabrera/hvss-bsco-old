import os
from inspect import stack
import polars as pl
import matplotlib.pyplot as plt
from tqdm import tqdm

###############################################################################

###############################################################################


def convert_headers(string):
    """
    Converts header string to list of column names.
    Searches for parts that start with "#" and have len > 3
        (normal format is "#%d:header"
    """
    split = string.split()
    headers = [s.split(":")[-1] for s in split if len(s) > 3 and s.startswith("#")]
    return headers


def cmc_model_to_params(model):
    """
    Converts cmc model name to tuple of params (N, rv, rg, Z).
    Drops "v2" if it is present.
    """
    split = model.strip(".tar.gz").split("_")
    v2_flag = 0
    if "v2" in split:
        v2_flag = 1
        split.remove("v2")
    N = float(split[0].strip("N")) / 1e5
    rv = float(split[1].strip("rv"))
    rg = float(split[2].strip("rg"))
    Z = float(split[3].strip("Z"))
    return N, rv, rg, Z, v2_flag


###############################################################################

CMC_CATALOG_PATH = "/hildafs/projects/phy200025p/share/catalog_files"

###############################################################################

# Iterate over all clusters in catalog directory
catalog_files = os.listdir(CMC_CATALOG_PATH)
count = 0
for model in tqdm(catalog_files):
    if model.endswith(".tar.gz"):
        count += 1
        print(model)
        # Try to load initial.bh.dat
        try:
            bh_path = "/".join((CMC_CATALOG_PATH, model, "initial.bh.dat"))
            bh_data = pl.read_csv(
                bh_path,
                has_header=False,
                sep=" ",
                skip_rows=1,
                dtypes=[
                    pl.datatypes.UInt32,
                    pl.datatypes.Float32,
                    pl.datatypes.UInt32,
                    pl.datatypes.UInt32,
                    pl.datatypes.UInt32,
                    pl.datatypes.UInt32,
                    pl.datatypes.UInt32,
                    pl.datatypes.UInt32,
                    pl.datatypes.UInt32,
                    pl.datatypes.UInt32,
                    pl.datatypes.UInt32,
                    pl.datatypes.UInt32,
                    pl.datatypes.Float32,
                ],
            )
            with open(bh_path) as file:
                bh_headers = file.readline()
            bh_data.columns = convert_headers(bh_headers)

            dyn_path = "/".join((CMC_CATALOG_PATH, model, "initial.dyn.dat"))
            dyn_data = pl.read_csv(
                dyn_path,
                has_header=False,
                columns=[7],
                new_columns=["r_c"],
                sep=" ",
                skip_rows=2,
            )
            #bh_data = bh_data.hstack(dyn_data)
        except Exception as e:
            print(
                "\t{file}, {line}: error in loading initial.bh.dat or initial.dyn.dat: {e}".format(
                    file=stack()[0][1].split("/")[-1], line=stack()[0][2], e=e
                )
            )
            print("\tcontinuing...")
            continue
        # Load intial.conv.sh (for converting TotalTime to Myr)
        # Copied from cmctools/cmctoolkit.py
        conv_path = "/".join((CMC_CATALOG_PATH, model, "initial.conv.sh"))
        f = open(conv_path, "r")
        convfile = f.read().split("\n")
        f.close()
        timeunitsmyr = float(convfile[19][13:])

        # The idea is that Nbh,tot for all clusters should be roughly split in time in two eras:
        # a. A rapid climb to Nbh_max after the first stars begin to collapse
        # b. An exponential decay to 0+ as bhs are ejected from the cluster.
        # A cluster is considered core-collapsed when Nbh,tot has decayed to < 10,
        # so t_cc is assigned by moving backwards through the timesteps until Nbh,tot >= 10 bhs.
        # There are some clusters that dissipate before any bhs are formed; these are definitely not core-collapsed.

        cci = -1
        if bh_data.max()[0, "Nbh,tot"] != 0:
            while bh_data[cci, "Nbh,tot"] < 10:
                cci -= 1
        cc_data_temp = bh_data[cci].clone()
        N, rv, rg, Z, _ = cmc_model_to_params(model)
        cc_data_temp = cc_data_temp.with_columns(
            [
                pl.lit(N).alias("N"),
                pl.lit(rv).alias("rv"),
                pl.lit(rg).alias("rg"),
                pl.lit(Z).alias("Z"),
                pl.when(cci == -1).then(pl.lit(0)).otherwise(pl.lit(1)).alias("ccflag"),
                (pl.col("TotalTime") * timeunitsmyr).alias("TotalTime")
            ]
        )
        cc_data_temp = cc_data_temp.select(
            [
                "N",
                "rv",
                "rg",
                "Z",
                "ccflag",
                "tcount",
                "TotalTime",
            ]
        )
        if count == 1:
            cc_data = cc_data_temp
        else:
            cc_data = cc_data.vstack(cc_data_temp)
print("{count} models processed.".format(count=count))
print("cc_data", cc_data)
cc_data.write_csv("cmc_core_collapsed.dat")
