"""Verify packaged release tables, stationarity scope, figures and file checksums."""
from pathlib import Path
import hashlib
import json
import pandas as pd

root = Path(__file__).resolve().parent
p = root / 'supplementary_data'
expected = [207, 2151, 1242, 621, 621, 5796, 621, 938]
for n, count in enumerate(expected, 1):
    files = list(p.glob(f'*S{n}_' + ('grafted_segments.csv' if n == 8 else '*.csv')))
    assert len(files) == 1, (n, files)
    assert len(pd.read_csv(files[0])) == count, (n, count)

stationarity = pd.read_csv(root / 'revision/equilibration_out/equilibration_per_replica.csv')
assert len(stationarity) == 621
assessed = stationarity[stationarity.status == 'ok']
assert len(assessed) == 618
flagged = assessed.iface_nonstationary.astype(str).str.lower().eq('true').sum()
assert flagged == 289, flagged

for entry in json.loads((root / 'figures/manifest.json').read_text()):
    for extension in ('png', 'pdf'):
        file = root / 'figures' / entry[extension]
        assert hashlib.sha256(file.read_bytes()).hexdigest() == entry[extension + '_sha256'], file
for entry in json.loads((root / 'SHA256SUMS.json').read_text()):
    file = root / entry['file']
    assert hashlib.sha256(file.read_bytes()).hexdigest() == entry['sha256'], file
print('Verified S1–S8 rows:', ' / '.join(map(str, expected)))
print('Verified interface stationarity: 289 / 618; release boundary: 621.')
print('Verified all packaged file and final-figure checksums; no trajectories read.')
