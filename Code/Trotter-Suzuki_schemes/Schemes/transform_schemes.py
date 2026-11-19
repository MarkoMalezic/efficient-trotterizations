"""Transform raw scheme files (schemes_q#.out) into other notations and formats.

Reads the scheme parameters a, b and the error/efficiency tables from the output of the
C++ program (see README-Trotter.md), and saves them together with
  - the "standard" parameters (notation of Omelyan et al., Eq. (8) of the paper),
  - the ramp parameters c, d (Eq. (9)),
  - the distance from the origin |(c, d) - 1/(2q)|,
as Numpy arrays and optionally in an HDF5 file. See README-Schemes.md for the file layout.

Usage:
    python transform_schemes.py path/to/schemes_q4.out                       # long double, saved next to the input
    python transform_schemes.py path/to/schemes_q4.out --precision double --out-dir some/dir/
    python transform_schemes.py --all                                        # regenerate every folder in Trotter-Suzuki/

Requires numpy, and h5py for the HDF5 output.
Quad precision is not supported, since numpy's long double is the 80-bit extended type on x86.
"""
import argparse
import glob
import os
import re

import numpy as np

PRECISIONS = {"double": np.float64, "long_double": np.longdouble}

_BLOCK = re.compile(r"(Scheme evaluation results:.*?)(?=Scheme evaluation results:|\Z)", re.DOTALL)
_CYCLES = re.compile(r"No\. cycles:\s*(\d+)")
_VECTORS = re.compile(r"Vector a:\s*([^\n]+)\s*\nVector b:\s*([^\n]+)")
_TABLE_ROW = re.compile(r"\|[ ]*(\d+)[ ]*\|[ ]*([0-9e\+\-\.]+|Infinity)[ ]*\|[ ]*([0-9e\+\-\.]+|Infinity)[ ]*\|")


def parse_schemes(filename):
  """Read all schemes in a schemes_q#.out file.

  Returns (q, a, b, err2, eff) with the parameters a, b and the squared errors and efficiencies
  per order (2, 4, 6, 8) as long double arrays with one row per scheme.
  An efficiency is set to infinity where the corresponding error vanishes exactly.
  """
  with open(filename, "r") as f:
    data = f.read()

  q = None
  a_vecs, b_vecs, err2s, effs = [], [], [], []
  for match in _BLOCK.finditer(data):
    block = match.group(1)
    vectors = _VECTORS.search(block)
    if not vectors:
      continue
    cycles = _CYCLES.search(block)
    if cycles:
      q = int(cycles.group(1))
    a_vecs.append([np.longdouble(x) for x in vectors.group(1).split(",")])
    b_vecs.append([np.longdouble(x) for x in vectors.group(2).split(",")])

    rows = _TABLE_ROW.findall(block)
    if not rows:
      continue
    err2 = [np.longdouble(row[1]) for row in rows]
    eff = [np.inf if e == 0.0 else np.longdouble(row[2]) for e, row in zip(err2, rows)]
    err2s.append(np.array(err2, dtype=np.longdouble))
    effs.append(np.array(eff, dtype=np.longdouble))

  as_array = lambda x: np.array(x, dtype=np.longdouble)
  return q, as_array(a_vecs), as_array(b_vecs), as_array(err2s), as_array(effs)


def to_standard(a, b, q):
  """Parameters of our implementation (one row per scheme) -> standard notation.

  Our parameters are stored in reverse order and multiplied by 1/2, except for the central one,
  which is a[0] for even q and b[0] for odd q.
  """
  a_std = a[:, ::-1] / 2.0
  b_std = b[:, ::-1] / 2.0
  if q % 2 == 0:
    a_std[:, q // 2] = a[:, 0]
  else:
    b_std[:, (q - 1) // 2] = b[:, 0]
  return a_std, b_std


def to_ramp(a_std, b_std, q):
  """Standard parameters -> ramp parameters c, d (q of each per scheme).

  Walks through the symmetric scheme from the outside in: c_i fills up the remaining part of the
  A step a_i, d_i the remaining part of the B step b_i, and the second half mirrors the first.
  """
  n = a_std.shape[0]
  c = np.zeros((n, q), dtype=a_std.dtype)
  d = np.zeros((n, q), dtype=b_std.dtype)
  half = (q + 1) // 2
  c[:, 0] = a_std[:, 0]
  d[:, 0] = b_std[:, 0] - c[:, 0]
  for i in range(1, q):
    if i < half:
      ia, ib = i, i
    else:
      back = i - half
      ia, ib = q // 2 - back, q // 2 - 1 - back
    c[:, i] = a_std[:, ia] - d[:, i - 1]
    d[:, i] = b_std[:, ib] - c[:, i]
  return c, d


def transform(filename, precision="long_double"):
  """Parse a schemes_q#.out file and compute all notations. Returns (q, dict of arrays)."""
  q, a, b, err2, eff = parse_schemes(filename)
  if q is None:
    raise ValueError(f"No schemes found in {filename}")
  dtype = PRECISIONS[precision]
  a, b, err2, eff = (x.astype(dtype) for x in (a, b, err2, eff))

  a_std, b_std = to_standard(a, b, q)
  c, d = to_ramp(a_std, b_std, q)
  distances = np.linalg.norm(np.hstack((c, d)) - 1 / (2 * q), axis=1)

  return q, {
    "parameters_a": a,
    "parameters_b": b,
    "errors_per_order": err2,
    "efficiencies_per_order": eff,
    "standard_parameters_a": a_std,
    "standard_parameters_b": b_std,
    "parameters_c": c,
    "parameters_d": d,
    "distances": distances,
  }


def save(filename, out_dir, precision="long_double", numpy=True, hdf5=True):
  """Transform a schemes_q#.out file and save the arrays as <name>_q#.npy and schemes_q#.h5 in out_dir."""
  q, arrays = transform(filename, precision)
  os.makedirs(out_dir, exist_ok=True)
  if numpy:
    for name, arr in arrays.items():
      np.save(os.path.join(out_dir, f"{name}_q{q}.npy"), arr)
  if hdf5:
    import h5py
    with h5py.File(os.path.join(out_dir, f"schemes_q{q}.h5"), "w") as hf:
      for name, arr in arrays.items():
        hf.create_dataset(name, data=arr, dtype=arr.dtype)
      hf.attrs["number of cycles"] = q
      hf.attrs["dtype"] = str(arrays["parameters_a"].dtype)
  return q, arrays


def all_jobs(root):
  """(input .out file, output directory, precision) for every scheme folder below root (= Trotter-Suzuki/)."""
  jobs = []
  # All_schemes/<precision>/order_n#/cycles_q#/schemes_q#.out -> the same folder
  for precision in PRECISIONS:
    for f in sorted(glob.glob(os.path.join(root, "All_schemes", precision, "order_n*", "cycles_q*", "schemes_q*.out"))):
      jobs.append((f, os.path.dirname(f), precision))
  # Recommended_schemes/order_n#/cycles_q#/schemes_q#.out -> a subfolder per precision
  for f in sorted(glob.glob(os.path.join(root, "Recommended_schemes", "order_n*", "cycles_q*", "schemes_q*.out"))):
    for precision in PRECISIONS:
      jobs.append((f, os.path.join(os.path.dirname(f), precision), precision))
  return jobs


if __name__ == "__main__":
  ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
  ap.add_argument("files", nargs="*", help="schemes_q#.out files to transform")
  ap.add_argument("--precision", default="long_double", choices=list(PRECISIONS))
  ap.add_argument("--out-dir", help="output directory (default: the directory of each input file)")
  ap.add_argument("--no-hdf5", action="store_true", help="only save the Numpy arrays")
  ap.add_argument("--all", action="store_true",
                  help="regenerate all folders in Trotter-Suzuki/ (All_schemes and Recommended_schemes, double and long double)")
  args = ap.parse_args()

  if args.all:
    jobs = all_jobs(os.path.join(os.path.dirname(os.path.abspath(__file__)), "Trotter-Suzuki"))
  else:
    if not args.files:
      ap.error("give at least one schemes_q#.out file, or --all")
    jobs = [(f, args.out_dir or os.path.dirname(os.path.abspath(f)), args.precision) for f in args.files]

  for f, out_dir, precision in jobs:
    q, arrays = save(f, out_dir, precision, hdf5=not args.no_hdf5)
    print(f"q = {q:2d}, {len(arrays['parameters_a']):5d} schemes, {precision:11s} -> {out_dir}")
