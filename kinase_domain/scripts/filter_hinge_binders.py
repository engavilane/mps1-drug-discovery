"""Filter PLIP results to hinge binders only."""
import pandas as pd
from pathlib import Path
import sys

input_csv  = sys.argv[1]
output_csv = sys.argv[2]

df     = pd.read_csv(input_csv)
hinge  = df[df['is_hinge_binder'] == True]
print(f"Hinge binders: {len(hinge)}/{len(df)}")
Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
hinge.to_csv(output_csv, index=False)
print(f"Saved → {output_csv}")
