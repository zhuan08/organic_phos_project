import os
import ase.io
from ase import Atoms
from ase.optimize import BFGS
from tblite.ase import TBLite
from ase.calculators.gamess_us import GAMESSUS
from ase.calculators.psi4 import Psi4
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDistGeom
import pandas as pd
IPythonConsole.ipython_3d = True

CHARGE = 0
METHOD = 'wb97m-v'
BASIS = 'def2-tzvppd'
additional_options = {
        'scf_type': 'df',
        'df_basis_scf': 'def2-universal-jkfit',
        'e_convergence': 1e-8,
        'd_convergence': 1e-8,
        'dft_radial_points': 99,
        'dft_spherical_points': 590,
        'dft_pruning_scheme': 'robust',
    }
nthreads = 1
memory = '90 GB'
psi4_calc_singlet = Psi4(method=METHOD, basis=BASIS, charge=CHARGE,
        multiplicity=1, reference = 'uks',
        num_threads = nthreads, memory = memory)
psi4_calc_singlet.psi4.set_options(additional_options)
psi4_calc_triplet = Psi4(method=METHOD, basis=BASIS, charge=CHARGE,
                  multiplicity=3, reference='uks',
                  num_threads = nthreads, memory = memory)
psi4_calc_triplet.psi4.set_options(additional_options)
tblite_calc_singlet = TBLite(multiplicity=1)
tblite_calc_triplet = TBLite(multiplicity=3)
gamess_calc_singlet = GAMESSUS(contrl={'mult': 1}, label='molecule',
                    command='rungms PREFIX.inp 30Jun2020R1 > PREFIX.log 2> PREFIX.err')
gamess_calc_triplet = GAMESSUS(contrl={'mult': 3}, label='molecule',
                    command='rungms PREFIX.inp 30Jun2020R1 > PREFIX.log 2> PREFIX.err')

geom_dir_name = 'new_organic_phos_geometries'
try:
    os.mkdir(geom_dir_name)
except FileExistsError:
    pass

test_mol_id = 'abp'
test_smile = "Nc1ccccc1C(=O)c1ccccc1"

def optimize_geometry(mol_id, smile, calc_triplet):
    geom_path = os.path.join(geom_dir_name, f'{mol_id}.xyz')
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
        atom.calc = calc_triplet
        opt = BFGS(atom, logfile=None, trajectory=None)
        opt.run(fmax=0.05)
        ase.io.write(filename=geom_path, images=atom)

def st_gap_calculate(path, calc_singlet, calc_triplet):
    atom = ase.io.read(filename=path)
    diff_energy = 0
    for multiplicity in ['singlet', 'triplet']:
        if multiplicity == 'singlet':
            atom.calc = calc_singlet
            energy = atom.get_potential_energy()
            diff_energy -= energy
        if multiplicity == 'triplet':
            atom.calc = calc_triplet
            energy = atom.get_potential_energy()
            diff_energy += energy
    return(diff_energy)

# Test using abp.xyz
if __name__ == "__main__":
    optimize_geometry(mol_id=test_mol_id, smile=test_smile, calc_triplet=tblite_calc_triplet)
    st_gap = st_gap_calculate(path=f"new_organic_phos_geometries/{test_mol_id}.xyz", calc_singlet=tblite_calc_singlet, calc_triplet=tblite_calc_triplet)
    print(st_gap)
