from pathlib import Path
from vina import Vina

v = Vina (sf_name='vina')

receptor_path = 'data/receptor/receptor.pdbqt'
ligands_dir = Path('data/ligands/pdbqt')

v.set_receptor('receptor.pdbqt')
v.set_ligand_from_file('ligand.pdbqt')
v.compute_vina_maps(center=[ -34.48, -15.66, -10.38], box_size=[20, 33, 21]) # Coordinates obtained through PyMOL, by using command line "centerofmass native_ligand.pdb" and by visual grid boxing
v.dock(exhaustiveness=16, n_poses=5)

results_dir = Path ('docking/results')
results_dir.mkdir(parents=True, exist_ok=True)
v.write_poses(str(results_dir / f'{ligand_name}_out.pdbqt', n_poses=5)
