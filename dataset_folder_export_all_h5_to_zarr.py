# dataset_folder_export_all_h5_to_zarr.py
# Converts a BDV-style fused H5 dataset into per-timepoint/view Zarr arrays.
# Output example: t00000_s00.zarr

import sys
import os
import re
import h5py
import numpy as np
import zarr
from numcodecs import Blosc

def parse_chunks(s, ndim=3):
    try:
        vals = tuple(int(x.strip()) for x in s.split(","))
        if len(vals) != ndim:
            raise ValueError
        return vals
    except Exception:
        return (64, 256, 256)

def get_dtype(outbits):
    return np.uint8 if str(outbits).startswith("8") else np.uint16

outbits = "16"
chunks = (64, 256, 256)

if len(sys.argv) > 1:
    if os.path.exists(sys.argv[1]):
        os.chdir(sys.argv[1])

if len(sys.argv) > 2:
    outbits = sys.argv[2]

if len(sys.argv) > 3:
    chunks = parse_chunks(sys.argv[3])

dtype = get_dtype(outbits)
compressor = Blosc(cname="zstd", clevel=5, shuffle=Blosc.BITSHUFFLE)

for root, dirs, files in os.walk("."):
    for filename in sorted(files):
        if not re.search(r"\.h5$", filename):
            continue

        print("Exporting Zarrs from file:", filename)

        with h5py.File(filename, "r") as f:
            for this_group_name in list(f.keys()):
                if not this_group_name.startswith("t0"):
                    continue

                print(" Working on group:", this_group_name)
                subgroup_names = list(f[f"/{this_group_name}/"].keys())

                for this_name in subgroup_names:
                    print("  subgroup:", this_name)

                    in_path = f"/{this_group_name}/{this_name}/0/cells"
                    if in_path not in f:
                        print("   skipping missing dataset:", in_path)
                        continue

                    out_name = f"{this_group_name}_{this_name}.zarr"

                    if os.path.exists(out_name):
                        print("   skipping existing:", out_name)
                        continue

                    arr = np.asarray(f[in_path], dtype=dtype)

                    z = zarr.open(
                        out_name,
                        mode="w",
                        shape=arr.shape,
                        dtype=arr.dtype,
                        chunks=tuple(min(c, s) for c, s in zip(chunks, arr.shape)),
                        compressor=compressor,
                    )
                    z[:] = arr
                    z.attrs["source_file"] = filename
                    z.attrs["source_dataset"] = in_path
                    z.attrs["axes"] = ["z", "y", "x"]
                    z.attrs["bit_depth"] = int(outbits)

                    print("   wrote:", out_name, arr.shape, arr.dtype)
