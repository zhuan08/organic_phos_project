import os
import ase.io
from ase import Atoms
from ase.optimize import BFGS
from tblite.ase import TBLite
from rdkit import Chem
from rdkit.Chem import rdDistGeom

def optimize_geometry(geom_path, mol_id, smile, triplet_setup):
    # First, check if there's already a geometry saved, and if so, just load it
    if os.path.exists(geom_path):
        atom = ase.io.read(filename=geom_path)
    else:
        molecule = Chem.MolFromSmiles(smile)
        if molecule is None:
            raise ValueError(f'MolFromSmiles returned None on {mol_id}')
        molecule = Chem.AddHs(molecule)
        rdDistGeom.EmbedMolecule(molecule)
        conf_mol = molecule.GetConformer()
        pos_mol = conf_mol.GetPositions()
        atom_mol = molecule.GetAtoms()
        atom_sym = []
        for char in atom_mol:
            atom_sym.append(char.GetSymbol())
        atom = Atoms(atom_sym, positions=pos_mol)
        triplet_setup(atom)
        opt = BFGS(atom, logfile=f'{mol_id}.opt', trajectory=f'{mol_id}.traj')
        opt.run(fmax=0.02)
        ase.io.write(filename=geom_path, images=atom)
    return atom

def st_gap_calculate(atom, singlet_setup, triplet_setup):
    diff_energy = 0
    for multiplicity in ['singlet', 'triplet']:
        if multiplicity == 'singlet':
            singlet_setup(atom)
            energy = atom.get_potential_energy()
            diff_energy -= energy
        if multiplicity == 'triplet':
            triplet_setup(atom)
            energy = atom.get_potential_energy()
            diff_energy += energy
    return diff_energy

def tblite_singlet_setup(atom):
    atom.calc = TBLite(multiplicity = 1)

def tblite_triplet_setup(atom):
    atom.calc = TBLite(multiplicity = 3)
# Test using abp.xyz
if __name__ == "__main__":
    geom_dir_name = 'new_organic_phos_geometries'
    try:
        os.mkdir(geom_dir_name)
    except FileExistsError:
        pass

    test_mol_id = 'abp'
    test_smile = "Nc1ccccc1C(=O)c1ccccc1"

    geom_path = os.path.join(geom_dir_name, f'tblite_{test_mol_id}.xyz')
    atom = optimize_geometry(geom_path=geom_path, mol_id=test_mol_id, smile=test_smile, triplet_setup=tblite_triplet_setup)
    st_gap = st_gap_calculate(atom, singlet_setup=tblite_singlet_setup, triplet_setup=tblite_triplet_setup)
    print(st_gap)
