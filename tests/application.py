#!/usr/bin/env python3
"""Exercise every registered problem and the fixed-dimension/module CLI contract."""
import csv
import json
import pathlib
import subprocess
import sys
import tempfile

exe, source, dimension = sys.argv[1:4]
dimension = int(dimension)
hydro, radiation, gravity = (value in ('1', 'ON', 'TRUE') for value in sys.argv[4:7])
runtime = sys.argv[7:]
problems = {
    'sod': ((1, 2, 3), hydro, 'HYDRO'),
    'kelvin-helmholtz': ((2, 3), hydro, 'HYDRO'),
    'rayleigh-taylor': ((3,), hydro, 'HYDRO'),
    'streaming': ((1, 2, 3), radiation, 'RADIATION'),
    'radiation-pulse': ((1, 2, 3), radiation, 'RADIATION'),
    'gravity-sphere': ((3,), gravity, 'GRAVITY'),
    'gravity-gaussian': ((3,), gravity, 'GRAVITY'),
    'collapse': ((3,), hydro and gravity, 'HYDRO' if not hydro else 'GRAVITY'),
    'polytrope': ((3,), hydro and gravity, 'HYDRO' if not hydro else 'GRAVITY'),
}

def run(*args, success=True, contains=None):
    result = subprocess.run([exe, *args, *runtime], text=True, capture_output=True, timeout=90)
    assert result.returncode == (0 if success else 1), (args, result.returncode, result.stdout, result.stderr)
    if contains:
        assert contains in result.stdout + result.stderr, (args, result.stdout, result.stderr)
    return result

run(success=False, contains='No problem specified')
run('--problem.name=', success=False)
run('--help', contains=f'octoII-{dimension}d')
for option in ('--mesh.ndim=3', '--ndim=3', '--problem.name=unknown'):
    run(option, success=False)
with tempfile.TemporaryDirectory(prefix='octoII-application-') as directory:
    root = pathlib.Path(directory)
    for name, (dimensions, enabled, module) in problems.items():
        if not enabled or dimension not in dimensions:
            expected = f'OCTOII_WITH_{module}=ON' if not enabled else 'supported dimensions:'
            run(f'--problem.name={name}', success=False, contains=expected)
            continue
        inputs = pathlib.Path(source) / 'examples' / (name + '.ini')
        if name == 'rayleigh-taylor':
            inputs = pathlib.Path(source) / 'examples/rayleigh-taylor-amr.ini'
        output = root / name
        run(f'--problem.name={name}', f'--config={inputs}', '--amr.enabled=off', '--amr.minLevel=0', '--mesh.cells=4', '--mesh.level=0',
            '--runtime.stopTime=' + ('0' if name.startswith('gravity-') else '1e-14'), '--output.every=1', '--output.enabled=on',
            f'--output.directory={output}', contains='Completed')
        rows = list(csv.DictReader((output / 'conservation.csv').open()))
        assert rows and rows[0]['step'] == '0', name
        if name in ('collapse', 'polytrope'):
            assert 'gas_gravity_energy_erg_grid' in rows[0], name
        reports = list(output.rglob('*.json'))
        assert reports, name
        for report in reports:
            data = json.loads(report.read_text())
            assert data['problem'] == name and data['ndim'] == dimension, data
    if hydro and radiation:
        # Resolve the final problem first, then apply that problem's defaults.
        first = root / 'first.ini'
        second = root / 'second.ini'
        first.write_text('[problem]\nname=sod\n[mesh]\ncells=8\n')
        second.write_text('[problem]\nname=streaming\n')
        # The INI can choose the problem; a command-line name still wins.
        run(f'--config={first}', '--mesh.cells=4', '--runtime.stopTime=0', '--output.enabled=off',
            f'--output.directory={root / "ini-only"}', contains='sod')
        run('--problem.name=streaming', f'--config={first}', f'--config={second}', '--mesh.cells=4', '--runtime.stopTime=0',
            '--output.enabled=off', contains='streaming')
        run('--problem.name=sod', f'--config={second}', '--runtime.stopTime=0', '--output.enabled=off', contains='sod')
print(f'{dimension}D: all problem availability, runtime selection, output identity and dimension rejection checks passed')
