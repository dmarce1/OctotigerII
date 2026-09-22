"""Exercise the normal executable, strict options and physical output."""
import csv
import math
from pathlib import Path
import subprocess
import sys

executable, examples, destination = map(Path, sys.argv[1:4])
runtime_args = sys.argv[4:]

def execute(*args, success=True):
    result = subprocess.run([str(executable), *args, *runtime_args], capture_output=True, text=True, timeout=120)
    if (result.returncode == 0) != success:
        raise RuntimeError(result.stdout + result.stderr)
    return result

def read_frame(path):
    with path.open() as stream:
        return [{k: float(v) for k, v in row.items()} for row in csv.DictReader(stream)]

execute('--not.an.option=1', success=False)
execute('--mesh.cells=7', success=False)
execute('--mesh.ndim=nan', success=False)
execute('--problem.name=gravity-sphere', '--mesh.periodic=on', success=False)
execute('--problem.name=gravity-sphere', '--mesh.ndim=2', success=False)
execute('--problem.name=streaming', '--radiation.light_speed_ratio=0', success=False)
execute('--runtime.max_steps=1', '--runtime.stop_time=1', '--output.enabled=off', success=False)

# Exact streaming translation at a reduced speed; physical flux must still
# be c*E, not c_hat*E. A pointwise relative L1 check catches wrong direction,
# missing evolution, and accidental use of the physical transport speed.
for ratio in (1.0, 0.25):
    folder = destination / str(ratio)
    execute('--config=' + str(examples / 'streaming.ini'), '--mesh.cells=64',
            '--runtime.stop_time=0.2', '--radiation.light_speed_ratio=' + str(ratio),
            '--output.directory=' + str(folder), '--output.every=10000')
    initial = read_frame(folder / 'frame_000000.csv')
    final = read_frame(folder / 'frame_000001.csv')
    length, speed, time = 6e10, 2.99792458e10, 0.2
    error, norm = 0.0, 0.0
    for row in final:
        distance = row['x_cm'] - (-3e10 + 0.25 * length + ratio * speed * time)
        distance -= round(distance / length) * length
        exact = 1e-6 + math.exp(-0.5 * (distance / (0.08 * length)) ** 2)
        error += abs(row['radiationEnergy'] - exact)
        norm += abs(exact)
        assert math.isclose(row['time_s'], time, rel_tol=1e-14)
        assert math.isclose(row['radiationFluxX'], speed * row['radiationEnergy'], rel_tol=1e-11)
    assert error / norm < 0.025, error / norm
    initial_energy = sum(row['radiationEnergy'] * row['dx_cm'] for row in initial)
    final_energy = sum(row['radiationEnergy'] * row['dx_cm'] for row in final)
    assert math.isclose(initial_energy, final_energy, rel_tol=2e-12)
    print(f'streaming ratio={ratio}: relative L1={error/norm:.6e}')

# Gravity-only output must contain source density and physical gravity.
folder = destination / 'gravity'
execute('--config=' + str(examples / 'gravity-gaussian.ini'), '--output.directory=' + str(folder))
rows = read_frame(folder / 'frame_000000.csv')
assert len(rows) == 8**3
assert all(row['potential'] < 0 and row['density'] > 0 for row in rows)
assert all(row['x_cm'] * row['accelerationX'] < 0 for row in rows)
assert all(math.isfinite(value) for row in rows for value in row.values())
print('Application options, cgs output, streaming translation and gravity output passed')
