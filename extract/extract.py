import parmap
import multiprocessing as mp

import rustics

###############################################################################

CHUNKSIZE = 1e4
MP_NPROCS = 128

mp_input = [
    "/".join((rustics.PATH_TO_DATA, name, rustics.FILE_REALZ))
    for name in rustics.MODEL_NAMES
]

extractor = rustics.EjectionExtractor(
    save_suffix=rustics.SUFFIX_EJECTIONS, chunksize=CHUNKSIZE
)

if MP_NPROCS == 1:
    # Don't run in parallel
    parmap.map(extractor.extract_ejections, mp_input, pm_parallel=False, pm_pbar=True)

else:
    # Run in parallel
    pool = mp.Pool(MP_NPROCS)
    parmap.map(extractor.extract_ejections, mp_input, pm_pool=pool, pm_pbar=True)
