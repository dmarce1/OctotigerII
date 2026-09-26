#!/usr/bin/env python3
"""Recreate equilibrium.inc from the pinned original Octo-Tiger paper dataset."""
import argparse
import hashlib
from pathlib import Path
import struct

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("original_binary", type=Path)
args = parser.parse_args()
data = args.original_binary.read_bytes()
expected = "9958252bae60348e08d4d043bf33954f75c43286133a99f49977e58d3b72da08"
if hashlib.sha256(data).hexdigest() != expected:
    parser.error("Input does not match the documented original equilibrium SHA256")
nr, nz, omega = struct.unpack_from("<iid", data)
assert (nr, nz, omega) == (100, 100, 0.5155532816213834)
assert len(data) == 16 + 16 * (2 * nr) * (2 * nz)
with Path(__file__).with_name("equilibrium.inc").open("w") as out:
    out.write("// Original Octo-Tiger SCF equilibrium; see DATA.md for provenance and license.\n")
    out.write("// Positive R,z quadrant, 100 x 100 cell centers; interleaved rho, internal energy.\n")
    for i in range(nr):
        for k in range(nz):
            rho, energy = struct.unpack_from("<dd", data, 16 + 16 * ((i + nr) * 2 * nz + k + nz))
            out.write("{" + rho.hex() + ", " + energy.hex() + "},\n")
