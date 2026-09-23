"""Exercise one compiled problem/dimension, strict options, and Silo time series."""
from pathlib import Path
import json
import math
import subprocess
import sys

executable, inputs, destination = map(Path, sys.argv[1:4])
problem, ndim = sys.argv[4:6]
runtime_args = sys.argv[6:]


def execute(*args, success=True):
    result = subprocess.run([str(executable), *args, *runtime_args],
                            capture_output=True, text=True, timeout=120)
    if (result.returncode == 0) != success:
        raise RuntimeError(result.stdout + result.stderr)
    return result


execute('--not.an.option=1', success=False)
execute('--output.format=csv', success=False)
execute('--mesh.cells=7', success=False)
# Even matching values cannot masquerade as runtime configuration choices.
execute('--problem.name=' + problem, success=False)
execute('--mesh.ndim=' + ndim, success=False)
execute('--radiation.lightSpeedRatio=0', success=False)
static = problem in ('gravity-sphere', 'gravity-gaussian')
if not static:
    execute('--runtime.maxSteps=1', '--runtime.stopTime=1e6', '--output.enabled=off',
            '--mesh.cells=4', '--mesh.level=0', success=False)
if problem in ('gravity-sphere', 'gravity-gaussian', 'collapse'):
    execute('--mesh.periodic=on', '--mesh.cells=4', '--mesh.level=0',
            '--runtime.stopTime=0', '--output.enabled=off')

result = execute('--config=' + str(inputs), '--mesh.cells=4', '--mesh.level=0',
                 '--runtime.stopTime=' + ('0' if static else '0.001'),
                 '--output.directory=' + str(destination), '--output.every=10000')
assert 'Completed' in result.stdout and problem in result.stdout and ndim + 'D' in result.stdout
frames = (destination / 'frames.visit').read_text().splitlines()
assert frames == [f'frame_{i:06d}.silo' for i in range(len(frames))]
assert len(frames) == (1 if static else 2)
assert all((destination / frame).is_file() for frame in frames)
print('Compiled problem/dimension, application options and Silo time series passed')

report = json.loads((destination / 'analytic-errors.json').read_text())
has_reference = problem in ('sod', 'gravity-sphere', 'gravity-gaussian', 'streaming', 'collapse')
assert report['status'] == ('available' if has_reference else 'unavailable')
assert report['problem'] == problem and report['ndim'] == int(ndim)
assert report['sampling'] == 'cell-center'
assert report['schemaVersion'] == 3
if has_reference:
    assert report['cells'] == 4 ** int(ndim)
    assert report['fields']
    for field in report['fields']:
        for norm in ('L1', 'L2', 'Linf'):
            absolute = field['absolute' + norm]
            reference = field['reference' + norm]
            relative = field[norm]
            assert absolute >= 0 and reference >= 0
            if reference:
                assert math.isclose(relative, absolute / reference, rel_tol=1e-14)
            else:
                assert relative is None
    execute('--mesh.cells=4', '--mesh.level=0', '--runtime.stopTime=' + ('0' if static else '0.05'),
            '--output.enabled=off', '--verification.analytic=on', '--verification.relativeL1Tolerance=0', success=False)
else:
    assert report['reason'] and not report['fields']
    execute('--output.enabled=off', '--verification.analytic=on', success=False)
execute('--verification.analytic=invalid', success=False)
execute('--verification.relativeL1Tolerance=nan', success=False)
execute('--verification.analytic=off', '--verification.relativeL1Tolerance=0.1', success=False)
print('Analytic report, availability, and failing accuracy gates passed')

if problem in ('gravity-sphere', 'gravity-gaussian', 'collapse'):
    assert report['referenceKind'] == 'direct'
    assert report['randomSeed'] == 5489 and report['randomGenerator'] == 'mt19937_64'
    args = ('--mesh.cells=4', '--mesh.level=0', '--runtime.stopTime=0',
            '--verification.directSamples=1', '--randomSeed=17',
            '--output.directory=' + str(destination))
    one = execute(*args)
    assert 'WARNING:' in one.stdout and 'cannot be estimated' in one.stdout
    sample = json.loads((destination / 'analytic-errors.json').read_text())
    assert sample['cells'] == 1 and sample['totalCells'] == 64
    assert sample['randomSeed'] == 17 and len(sample['targetIndices']) == 1
    assert sample['sampling'] == 'random-cell-center-without-replacement'
    assert all(f['samplingWarning'] and f['relativeL1HalfWidth95'] is None for f in sample['fields'])
    execute(*args)
    assert sample == json.loads((destination / 'analytic-errors.json').read_text())
