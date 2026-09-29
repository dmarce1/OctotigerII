#!/usr/problem/env python3
"""Reproduce the precision/quadrature experiment without changing the EOS library.

C++17 compiler required; mpmath required only for orders other than 20.
Generated sources/binaries and CSV measurements live in --output-dir.
The Sommerfeld column is a leading low-temperature comparison, not an EOS
substitute; only interpret it when kT/EF is small.
"""
import argparse
from pathlib import Path
import re
import subprocess

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--output-dir', type=Path, default=Path('build-helmholtz/precision-experiment'))
p.add_argument('--orders', type=int, nargs='+', default=[20, 40])
p.add_argument('--compiler', default='c++')
a = p.parse_args()
root = Path(__file__).resolve().parents[1]
out = a.output_dir.resolve()
out.mkdir(parents=True, exist_ok=True)
if any(n < 20 or n % 2 for n in a.orders):
    p.error('Orders must be even and at least 20')

def source_variant(text, extended):
    namespace = 'extended' if extended else 'baseline'
    result = []
    for line in text.splitlines():
        if line.lstrip().startswith('//'):
            result.append(line)
            continue
        line = line.replace('octotigerII::helmholtz::detail', 'experiment::' + namespace)
        line = line.replace('"direct.hpp"', f'"direct_{namespace}.hpp"')
        line = line.replace('"constants.hpp"', f'"constants_{namespace}.hpp"')
        if extended:
            line = re.sub(r'\bdouble\b', 'long double', line)
            line = re.sub(r'(?<![\w.])(?:\d+\.\d*(?:[eE][+-]?\d+)?|\d+[eE][+-]?\d+)(?![\w.])',
                          lambda m: m[0] + 'L', line)
        result.append(line)
    text = '\n'.join(result) + '\n'
    if extended:
        # Promoting already-rounded double coefficients would not test full
        # long-double quadrature. Recover Timmes's original decimal literals.
        pat = r'(//\s*data\s+(xg|wg)\s*\(\s*(\d+)\s*\)\s*/\s*([^/]+)\s*/\s*\n\s*)(?:xg|wg)\[\d+\] = [^;]+;'
        text = re.sub(pat, lambda m: m[1] + f'{m[2]}[{m[3]}] = ' +
                      m[4].strip().replace('d', 'e') + 'L;', text)
        text = text.replace('3.141592653589793L', '3.1415926535897932384L')
        text = text.replace('0.5772156649015329L', '0.577215664901532861L')
    # Expose the independent energy derivative for this experiment only.
    text = text.replace('n, nd, nt, ndt;', 'n, nd, nt, ndt, cvEnergy;')
    text = text.replace('w.dxneferddt + w.dxnpferddt};', 'w.dxneferddt + w.dxnpferddt, w.deepdt};')
    return text

def higher_order(text, order):
    if order == 20:
        return text
    import mpmath as mp
    mp.mp.dps = 60
    for name, family, half in [('dqleg020', 'legendre', True), ('dqlag020', 'laguerre', False)]:
        start = text.rfind('void ' + name + '(')
        end = text.index('\n}', start) + 2
        part = text[start:end]
        nodes, weights = mp.gauss_quadrature(order, family)
        if half:
            nodes, weights = list(nodes)[order//2:], list(weights)[order//2:]
        else:
            weights = [w*mp.exp(x) for w, x in zip(weights, nodes)]
        count = len(nodes)
        part = re.sub(r'long double ([wx]g)\[\d+\]', rf'long double \1[{count+1}]', part)
        part = re.sub(r'^\s*[wx]g\[\d+\] = [^;]+;', '', part, flags=re.M)
        part = part.replace('j <= ' + ('10' if half else '20'), f'j <= {count}')
        init = f'\n// Experiment: {order}-point coefficients computed at 60 decimal digits.\n'
        for key, vals in [('xg', nodes), ('wg', weights)]:
            for i, value in enumerate(vals, 1):
                init += f'{key}[{i}] = {mp.nstr(value,40)}L;\n'
        anchor = '    center =' if half else '    result = 0.0L;'
        pos = part.index(anchor)
        part = part[:pos] + init + part[pos:]
        text = text[:start] + part + text[end:]
    return text

for extended in (False, True):
    ns = 'extended' if extended else 'baseline'
    for filename in ('direct.hpp', 'constants.hpp', 'direct.cpp'):
        target = out / filename.replace('.', '_' + ns + '.', 1)
        target.write_text(source_variant((root/filename).read_text(), extended))

probe = r'''
#include "direct_baseline.hpp"
#include "direct_extended.hpp"
#include "constants_extended.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
int main() {
    using namespace experiment::extended;
    std::cout << "D,T,double_bits,long_double_bits,double_cv_entropy,double_cv_energy,"
                 "long_double_cv_entropy,long_double_cv_energy,long_double_entropy,"
                 "sommerfeld_cv,kT_over_fermi_energy\n" << std::setprecision(18);
    for (long double d : {1e6L,1e9L,1e12L,1e15L}) {
        for (long double t : {1e3L,1e4L,1e6L,1e8L}) {
            auto b = experiment::baseline::direct(d,t);
            auto q = direct(d,t);
            const long double x=hbar*std::cbrt(3*pi*pi*avo*d)/(me*clight);
            const long double ef=me*clight*clight*x*x/(std::sqrt(1+x*x)+1);
            // Density of states at EF: g(EF)=3n*Etotal/(pF^2 c^2).
            // Cv=(pi^2/3)*kB^2*T*g(EF)/D to leading Sommerfeld order.
            const long double cv=pi*pi*avo*kerg*kerg*t*std::sqrt(1+x*x)/(me*clight*clight*x*x);
            std::cout << d << ',' << t << ',' << std::numeric_limits<double>::digits << ','
                << std::numeric_limits<long double>::digits << ',' << t*b.st << ',' << b.cvEnergy
                << ',' << t*q.st << ',' << q.cvEnergy << ',' << q.s << ',' << cv << ',' << kerg*t/ef << '\n';
        }
    }
}
'''
(out/'probe.cpp').write_text(probe)
base = (out/'direct_extended.cpp').read_text()
for order in a.orders:
    variant = out/f'direct_extended_{order}.cpp'
    variant.write_text(higher_order(base, order))
    exe = out/f'probe_{order}'
    subprocess.run([a.compiler, '-O2', '-std=c++17', '-I'+str(out),
                    str(out/'direct_baseline.cpp'), str(variant), str(out/'probe.cpp'), '-o', str(exe)], check=True)
    result = subprocess.run([str(exe)], check=True, text=True, capture_output=True)
    csv = out/f'precision_{order}.csv'
    csv.write_text(result.stdout)
    print(f'{order}-point quadrature: {csv}')
    lines = result.stdout.splitlines()
    print(lines[0])
    print(next(line for line in lines if line.startswith('1000000000,10000,')))
