#!/usr/bin/env python3
"""Configure, build and test the supported problem/dimension combinations.

The problem manifests supply the matrix. Each build has a separate directory;
CMake prefix paths and compiler choices are inherited from the environment.
"""
import argparse
import os
from pathlib import Path
import re
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backend", choices=("serial", "hpx"), default="serial")
    parser.add_argument("--problem", action="append", help="Limit to this problem (repeatable)")
    parser.add_argument("--ndim", type=int, choices=(1, 2, 3), help="Limit to this dimension")
    parser.add_argument("--build-root", type=Path)
    parser.add_argument("--jobs", "-j", type=int, default=min(4, os.cpu_count() or 1))
    parser.add_argument("--build-type", default="Release", choices=("Release", "Debug", "RelWithDebInfo"))
    parser.add_argument("--unit-only", action="store_true")
    parser.add_argument("--cmake-arg", action="append", default=[], help="Extra CMake argument, e.g. --cmake-arg=-DHPX_DIR=/path")
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    source = Path(__file__).resolve().parents[1]
    root = (args.build_root or source / ".test-build" / args.backend / args.build_type.lower()).resolve()
    combinations = []
    names = set()
    for manifest in sorted((source / "bin").glob("*/*/CMakeLists.txt")):
        match = re.search(r"octo_problem\(([^\s]+)\s+DIMENSIONS\s+([123\s]+)DEFAULT_DIMENSION", manifest.read_text())
        if not match:
            raise RuntimeError(f"Cannot read problem dimensions from {manifest}")
        name, dimensions = match.groups()
        names.add(name)
        if args.problem and name not in args.problem:
            continue
        for dimension in map(int, dimensions.split()):
            if args.ndim is None or dimension == args.ndim:
                combinations.append((name, dimension))
    if args.problem and set(args.problem) - names:
        parser.error("Unknown problem(s): " + ", ".join(sorted(set(args.problem) - names)))
    if not combinations:
        parser.error("No supported problem/dimension combinations selected")
    for problem, dimension in combinations:
        build = root / problem / f"{dimension}d"
        print(f"Testing {problem} {dimension}D ({args.backend})", flush=True)
        subprocess.run(["cmake", "-S", str(source), "-B", str(build), *args.cmake_arg,
                        f"-DCMAKE_BUILD_TYPE={args.build_type}", "-DOCTOTIGERII_BUILD_TESTS=ON",
                        f"-DOCTOTIGERII_WITH_HPX={'ON' if args.backend == 'hpx' else 'OFF'}",
                        f"-DOCTOTIGERII_PROBLEM={problem}", f"-DOCTOTIGERII_NDIM={dimension}"], check=True)
        subprocess.run(["cmake", "--build", str(build), "--parallel", str(args.jobs)], check=True)
        command = ["ctest", "--test-dir", str(build), "--output-on-failure", "--parallel", str(args.jobs),
                   "--output-junit", str(build / "test-results.xml")]
        if args.unit_only:
            command += ["-L", "^unit$"]
        subprocess.run(command, check=True)
    print(f"Passed {len(combinations)} problem/dimension builds.", flush=True)


if __name__ == "__main__":
    main()
