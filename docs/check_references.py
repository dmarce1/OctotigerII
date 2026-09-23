"""Reject duplicate bibliography anchors and broken source citation keys."""
from pathlib import Path
import re

root = Path(__file__).resolve().parents[1]
bibliography = (root / 'BIBLIOGRAPHY.md').read_text()
anchors = re.findall(r'\{#(ref_\w+)\}', bibliography)
if len(anchors) != len(set(anchors)):
    raise SystemExit('Duplicate bibliography anchor')
known = set(anchors)
used = set()
for folder in ('include', 'src', 'bin', 'docs'):
    for path in (root / folder).rglob('*'):
        if path.suffix not in ('.cpp', '.hpp', '.md', '.dox') or 'generated' in path.parts:
            continue
        for key in re.findall(r'[@\\]ref\s+(ref_\w+)', path.read_text()):
            if key not in known:
                raise SystemExit(f'{path.relative_to(root)}: unresolved citation {key}')
            used.add(key)
missing = known - used
if missing:
    raise SystemExit('Bibliography entries without a code/manual citation: ' + ', '.join(sorted(missing)))
print(f'Bibliography: {len(known)} records, all citation keys resolve')
