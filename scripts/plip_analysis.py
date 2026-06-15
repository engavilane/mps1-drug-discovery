"""
PLIP-based protein-ligand interaction analysis.

Usage:
  python scripts/plip_analysis.py \
    --scores  docking/results/phase1_docking_scores.csv \
    --results docking/results \
    --output  analysis/interactions_plip/phase1 \
    --old_interactions analysis/interactions/interaction_analysis.csv

  python scripts/plip_analysis.py \
    --scores  docking/phase2_results/docking_scores.csv \
    --results docking/phase2_results \
    --output  analysis/interactions_plip/phase2
"""

import argparse
from plip.structure.preparation import PDBComplex
import pandas as pd
import numpy as np
from pathlib import Path
import tempfile
import os
import warnings
warnings.filterwarnings('ignore')


# Argparse 
ap = argparse.ArgumentParser()
ap.add_argument("--scores",
    default="docking/results/phase1_docking_scores.csv")
ap.add_argument("--results",
    default="docking/results")
ap.add_argument("--output",
    default="analysis/interactions_plip/phase1")
ap.add_argument("--receptor",
    default="data/receptor/receptor_clean.pdb")
ap.add_argument("--old_interactions",
    default=None)
args = ap.parse_args()

RECEPTOR_PDB = args.receptor
RESULTS_DIR  = Path(args.results)
SCORES_CSV   = args.scores
OUTPUT_DIR   = Path(args.output)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
HINGE        = {603, 605}

print(f"PLIP Interaction Analysis")
print(f"  Scores:   {SCORES_CSV}")
print(f"  Results:  {RESULTS_DIR}")
print(f"  Output:   {OUTPUT_DIR}")
print(f"  Receptor: {RECEPTOR_PDB}\n")


# Helpers 
def pdbqt_to_pdb_lines(pdbqt_path):
    lines    = []
    atom_num = 1
    with open(pdbqt_path) as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            at = line[77:79].strip() if len(line) > 77 else ""
            if at.startswith("H"):
                continue
            try:
                name_a = line[12:16].strip()
                x      = float(line[30:38])
                y      = float(line[38:46])
                z      = float(line[46:54])
            except (ValueError, IndexError):
                continue
            lines.append(
                f"HETATM{atom_num:5d} {name_a:<4s} "
                f"LIG A{999:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}"
                f"  1.00  0.00          \n"
            )
            atom_num += 1
    return lines

def create_complex_pdb(receptor_pdb, ligand_lines, out_path):
    with open(out_path, "w") as out:
        with open(receptor_pdb) as rec:
            for line in rec:
                if line.startswith(("ATOM", "HETATM")):
                    out.write(line)
        for line in ligand_lines:
            out.write(line)
        out.write("END\n")

def run_plip_on_complex(complex_pdb):
    mol = PDBComplex()
    mol.load_pdb(complex_pdb)
    skip = {"HOH","EDO","GOL","SO4","MG","CL","PEG","ACT"}
    for lig in mol.ligands:
        if lig.hetid == "LIG" or lig.hetid not in skip:
            mol.characterize_complex(lig)
            key = f"{lig.hetid}:{lig.chain}:{lig.position}"
            if key in mol.interaction_sets:
                return mol.interaction_sets[key]
    return None

def count_hinge_hbonds(ix):
    if ix is None:
        return 0, 0, []
    pdon, ldon = 0, 0
    details    = []
    for hb in ix.hbonds_pdon:
        if hb.resnr in HINGE:
            pdon += 1
            details.append(
                f"{hb.restype}{hb.resnr}(pdon)"
                f"d={hb.distance_ad:.2f}"
                f"a={hb.angle:.1f}"
            )
    for hb in ix.hbonds_ldon:
        if hb.resnr in HINGE:
            ldon += 1
            details.append(
                f"{hb.restype}{hb.resnr}(ldon)"
                f"d={hb.distance_ad:.2f}"
                f"a={hb.angle:.1f}"
            )
    return pdon, ldon, details

def count_other_interactions(ix):
    if ix is None:
        return 0, 0, 0, 0, set()
    hydro, pi, salt, water = 0, 0, 0, 0
    res_hit = set()
    for hc in ix.hydrophobic_contacts:
        if hc.resnr in HINGE:
            hydro += 1
            res_hit.add(hc.resnr)
    for ps in ix.pistacking:
        if ps.resnr in HINGE:
            pi += 1
            res_hit.add(ps.resnr)
    for sb in (ix.saltbridge_lneg + ix.saltbridge_pneg):
        if sb.resnr in HINGE:
            salt += 1
            res_hit.add(sb.resnr)
    for wb in ix.water_bridges:
        if wb.resnr in HINGE:
            water += 1
            res_hit.add(wb.resnr)
    return hydro, pi, salt, water, res_hit


# Load scores 
print("Loading docking scores...")
scores_df = pd.read_csv(SCORES_CSV)
print(f"  {len(scores_df)} ligands to analyse\n")
print("=" * 65)

results = []
failed  = []


# Main loop 
for idx, row in scores_df.iterrows():
    ligand_name = row["ligand"]
    score       = row["affinity_kcal"]

    pdbqt_path = RESULTS_DIR / f"{ligand_name}_out.pdbqt"
    if not pdbqt_path.exists():
        pdbqt_path = RESULTS_DIR / f"{ligand_name}.pdbqt"
    if not pdbqt_path.exists():
        failed.append(ligand_name)
        print(f"  {ligand_name:<30} ✗ PDBQT not found")
        continue

    tmp_path = None
    try:
        ligand_lines = pdbqt_to_pdb_lines(str(pdbqt_path))
        if not ligand_lines:
            failed.append(ligand_name)
            print(f"  {ligand_name:<30} ✗ No ligand atoms")
            continue

        with tempfile.NamedTemporaryFile(
            suffix=".pdb", delete=False, mode="w"
        ) as tmp:
            tmp_path = tmp.name

        create_complex_pdb(
            RECEPTOR_PDB, ligand_lines, tmp_path
        )

        ix = run_plip_on_complex(tmp_path)

        pdon, ldon, details = count_hinge_hbonds(ix)
        hydro, pi, salt, water, res_hit = \
            count_other_interactions(ix)

        total    = pdon + ldon
        is_hinge = total > 0

        hinge_res = sorted(list(res_hit))
        for hb in ix.hbonds_pdon if ix else []:
            if hb.resnr in HINGE and hb.resnr not in hinge_res:
                hinge_res.append(hb.resnr)
        for hb in ix.hbonds_ldon if ix else []:
            if hb.resnr in HINGE and hb.resnr not in hinge_res:
                hinge_res.append(hb.resnr)
        hinge_res = sorted(set(hinge_res))

        results.append({
            "ligand":          ligand_name,
            "affinity_kcal":   score,
            "hbond_total":     total,
            "hbond_pdon":      pdon,
            "hbond_ldon":      ldon,
            "hydrophobic":     hydro,
            "pi_stack":        pi,
            "salt_bridge":     salt,
            "water_bridge":    water,
            "hinge_residues":  str(hinge_res),
            "hbond_details":   " | ".join(details),
            "is_hinge_binder": is_hinge,
        })

        tag = "✓ hinge" if is_hinge else "· no hinge"
        print(f"  {ligand_name:<30} "
              f"HB:{total} "
              f"Hyd:{hydro} "
              f"π:{pi}  "
              f"{tag}")

    except Exception as e:
        import traceback
        failed.append(ligand_name)
        print(f"  {ligand_name:<30} ✗ {e}")
        traceback.print_exc()

    finally:
        if tmp_path and os.path.exists(tmp_path):
            os.unlink(tmp_path)


# Build DataFrame 
if not results:
    print("\n✗ No results — check errors above")
    exit(1)

df = pd.DataFrame(results)
df.to_csv(OUTPUT_DIR / "plip_interactions.csv", index=False)


# Summary 
hinge_n    = int(df["is_hinge_binder"].sum())
no_hinge_n = int((~df["is_hinge_binder"]).sum())

print("\n" + "=" * 65)
print("PLIP ANALYSIS COMPLETE")
print("=" * 65)
print(f"  Hinge binders (≥1 H-bond): {hinge_n}/{len(df)}")
print(f"  No hinge H-bond:           {no_hinge_n}/{len(df)}")
print(f"  Failed:                    "
      f"{len(failed)}/{len(scores_df)}")

print(f"\nTop 10 by affinity:")
print(f"{'Ligand':<30} {'Score':>7} "
      f"{'HBonds':>7} {'Hydro':>6} {'π':>4}")
print("-" * 55)
for _, r in df.nsmallest(10, "affinity_kcal").iterrows():
    print(f"  {r['ligand']:<28} "
          f"{r['affinity_kcal']:>7.3f} "
          f"{r['hbond_total']:>7} "
          f"{r['hydrophobic']:>6} "
          f"{r['pi_stack']:>4}")


# Compare with distance-based 
if args.old_interactions and \
        Path(args.old_interactions).exists():
    print(f"\nComparison with distance-based analysis:")
    old = pd.read_csv(args.old_interactions)

    if "gly605_contacts" in old.columns:
        old["old_hinge"] = (
            old["gly605_contacts"] +
            old["glu603_contacts"]
        ) > 0
    elif "is_hinge_binder" in old.columns:
        old["old_hinge"] = old["is_hinge_binder"]
    else:
        old["old_hinge"] = True

    merged = df.merge(
        old[["ligand", "old_hinge"]],
        on="ligand", how="inner"
    )

    agree  = (merged["is_hinge_binder"] ==
              merged["old_hinge"]).sum()
    gained = merged[
        merged["is_hinge_binder"] & ~merged["old_hinge"]
    ]["ligand"].tolist()
    lost   = merged[
        ~merged["is_hinge_binder"] & merged["old_hinge"]
    ]["ligand"].tolist()

    print(f"  Agreement:           {agree}/{len(merged)}")
    print(f"  Newly hinge binders: {len(gained)}")
    for l in gained:
        print(f"    + {l}")
    print(f"  Lost hinge binders:  {len(lost)}")
    for l in lost:
        print(f"    - {l}")

print(f"\nResults → {OUTPUT_DIR}/plip_interactions.csv")
