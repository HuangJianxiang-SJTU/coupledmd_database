from pathlib import Path
import pandas as pd
p=Path(__file__).resolve().parent/'supplementary_data'
expected=[207,2151,1242,621,621,5796,621,938]
for n,count in enumerate(expected,1):
 f=next(p.glob(f'*S{n}_'+('grafted_segments.csv' if n==8 else '*.csv')))
 assert len(pd.read_csv(f))==count,(n,count)
print('S1–S8 cohort boundary verified; no trajectories read.')
