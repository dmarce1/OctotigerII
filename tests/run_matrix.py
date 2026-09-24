#!/usr/bin/env python3
"""Build/test all dimensions, optionally across all eight module combinations."""
import argparse
import itertools
import os
from pathlib import Path
import subprocess

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--backend', choices=('serial', 'hpx'), default='serial')
p.add_argument('--build-root', type=Path)
p.add_argument('--jobs', '-j', type=int, default=min(4, os.cpu_count() or 1))
p.add_argument('--build-type', choices=('Release', 'Debug', 'RelWithDebInfo'), default='Release')
p.add_argument('--unit-only', action='store_true')
p.add_argument('--module-matrix', action='store_true')
p.add_argument('--cmake-arg', action='append', default=[])
a = p.parse_args()
if a.jobs < 1:
    p.error('--jobs must be positive')
source = Path(__file__).resolve().parents[1]
root = (a.build_root or source / '.test-build' / a.backend / a.build_type.lower()).resolve()
combinations = list(itertools.product(('ON', 'OFF'), repeat=3)) if a.module_matrix else [('ON',) * 3]
for flags in combinations:
    build = root / '-'.join(flags)
    subprocess.run(['cmake', '-S', str(source), '-B', str(build), *a.cmake_arg,
                    f'-DCMAKE_BUILD_TYPE={a.build_type}', '-DOCTOTIGERII_BUILD_TESTS=ON',
                    f'-DOCTOTIGERII_WITH_HPX={"ON" if a.backend == "hpx" else "OFF"}',
                    *(f'-DOCTOII_WITH_{m}={v}' for m, v in zip(('HYDRO', 'RADIATION', 'GRAVITY'), flags))], check=True)
    subprocess.run(['cmake', '--build', str(build), '-j', str(a.jobs)], check=True)
    command = ['ctest', '--test-dir', str(build), '--output-on-failure', '-j', str(a.jobs),
               '--output-junit', str(build / 'test-results.xml')]
    if a.unit_only:
        command += ['-L', '^unit$']
    subprocess.run(command, check=True)
print(f'Passed {len(combinations)} configurations, three dimensions each.')
