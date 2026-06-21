"""
PLIP interaction analysis for NTD TPR domain docking results.

Checks interactions with Hec1/Ndc80 interface residues
identified by Screpanti et al. 2011:
  104, 107, 108, 111 (TPR1)
  140, 141, 144, 145 (TPR2)
  173, 174, 177, 178 (TPR3)
"""

from plip.structure.preparation import PDBComplex
import pandas as pd
import numpy as np
from pathlib import Path
import tempfile
import os
import argparse
import warnings
warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser()
parser.add_argument("--scores",
    default="ntd_domain/docking/top20_results/top20_scores.csv")
parser.add_argument("--results",
    default="ntd_domain/docking/top20_results")
parser.add_argument("--output",
    default="ntd_domain/analysis/plip_top20")
parser.add_argument("--receptor",
    default="ntd_domain/data/receptor/4H7Y_receptor.pdbqt")
args = parser.parse_args()

RESULTS_DIR = Path(args.results)
OUTPUT_DIR  = Path(args.output)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Hec1 interface residues from Screpanti et al. 2011
HEC1_INTERFACE = {104, 107, 108, 111,
                  140, 141, 144, 145,
                  173, 174, 177, 178}

print("NTD PLIP Interaction Analysis")
print("=" * 55)
print(f"Interface residues: {sorted(HEC1_INTERFACE)}")

def pdbqt_to_pdb_lines(pdbqt_path):
    lines    = []
    atom_num = 1
    with open(pdbqt_path) as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
            at = line[77:79].strip() if len(line) > 77 else ''
            if at.startswith('H'):
                continue
            try:
                name_a = line[12:16].strip()
                x      = float(line[30:38])
                y      = float(line[38:46])
                z      = float(line[46:54])
            except:
                continue
            lines.append(
                f'HETATM{atom_num:5d} {name_a:<4s} '
                f'LIG A{999:4d}    '
                f'{x:8.3f}{y:8.3f}{z:8.3f}'
                f'  1.00  0.00          \n'
            )
            atom_num += 1
    return lines

def receptor_to_pdb_lines(pdbqt_path):
    lines = []
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                lines.append(f'{line[:66].rstrip():<66}\n')
    return lines

scores_df = pd.read_csv(args.scores)
print(f"Ligands to analyse: {len(scores_df)}")
print("=" * 55)

rec_lines = receptor_to_pdb_lines(args.receptor)
results   = []

for _, row in scores_df.iterrows():
    ligand_name = row['ligand']
    score       = row['affinity_kcal']

    pdbqt_path = RESULTS_DIR / f"{ligand_name}_out.pdbqt"
    if not pdbqt_path.exists():
        print(f"  {ligand_name}: PDBQT not found")
        continue

    lig_lines = pdbqt_to_pdb_lines(str(pdbqt_path))
    if not lig_lines:
        print(f"  {ligand_name}: no atoms extracted")
        continue

    with tempfile.NamedTemporaryFile(
        suffix='.pdb', delete=False, mode='w'
    ) as tmp:
        tmp_path = tmp.name
        for line in rec_lines:
            tmp.write(line)
        for line in lig_lines:
            tmp.write(line)
        tmp.write('END\n')

    try:
        mol = PDBComplex()
        mol.load_pdb(tmp_path)
        mol.analyze()

        hbonds      = 0
        hydrophobic = 0
        pi_stack    = 0
        interface_contacts = set()
        hbond_details = []

        for lig_key, interactions in mol.interaction_sets.items():
            for hb in interactions.hbonds_pdon:
                res_id = hb.resnr
                if res_id in HEC1_INTERFACE:
                    hbonds += 1
                    interface_contacts.add(res_id)
                    hbond_details.append(
                        f"{hb.restype}{res_id}"
                        f"(pdon)d={hb.distance_ad:.2f}"
                        f"a={hb.angle:.1f}"
                    )
            for hb in interactions.hbonds_ldon:
                res_id = hb.resnr
                if res_id in HEC1_INTERFACE:
                    hbonds += 1
                    interface_contacts.add(res_id)
                    hbond_details.append(
                        f"{hb.restype}{res_id}"
                        f"(ldon)d={hb.distance_ad:.2f}"
                        f"a={hb.angle:.1f}"
                    )
            for hc in interactions.hydrophobic_contacts:
                if hc.resnr in HEC1_INTERFACE:
                    hydrophobic += 1
                    interface_contacts.add(hc.resnr)
            for ps in interactions.pistacking:
                if ps.resnr in HEC1_INTERFACE:
                    pi_stack += 1
                    interface_contacts.add(ps.resnr)

        is_interface_binder = len(interface_contacts) > 0
        status = '✓ interface' if is_interface_binder else '· no contact'

        print(f"  {ligand_name:<25} "
              f"HB:{hbonds} "
              f"Hyd:{hydrophobic} "
              f"π:{pi_stack}  "
              f"{status}")

        results.append({
            'ligand':              ligand_name,
            'affinity_kcal':       score,
            'hbond_total':         hbonds,
            'hydrophobic':         hydrophobic,
            'pi_stack':            pi_stack,
            'interface_residues':  str(sorted(interface_contacts)),
            'hbond_details':       ' | '.join(hbond_details),
            'is_interface_binder': is_interface_binder,
        })

    except Exception as e:
        print(f"  {ligand_name}: PLIP error — {e}")
    finally:
        os.unlink(tmp_path)

results_df = pd.DataFrame(results)
results_df.to_csv(
    OUTPUT_DIR / "ntd_plip_interactions.csv", index=False
)

interface_binders = results_df[
    results_df['is_interface_binder'] == True
]

print(f"\n{'=' * 55}")
print(f"NTD PLIP ANALYSIS COMPLETE")
print(f"{'=' * 55}")
print(f"  Interface binders: "
      f"{len(interface_binders)}/{len(results_df)}")

print(f"\nTop candidates by affinity + interface contact:")
print(f"{'Ligand':<25} {'Score':>7} {'HB':>3} "
      f"{'Hyd':>4} {'Residues':>20}")
print("-" * 65)

for _, r in results_df.sort_values(
    'affinity_kcal'
).head(15).iterrows():
    marker = '✓' if r['is_interface_binder'] else ' '
    print(f"  {marker} {r['ligand']:<23} "
          f"{r['affinity_kcal']:>7.3f} "
          f"{r['hbond_total']:>3} "
          f"{r['hydrophobic']:>4} "
          f"{r['interface_residues']:>20}")

print(f"\nResults → {OUTPUT_DIR}/ntd_plip_interactions.csv")
